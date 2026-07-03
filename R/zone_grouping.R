#' Validate and Normalize External Zone Assignments
#'
#' @param zone_assignments A data.frame with `tract_id`, `zone_id`, and
#'   optionally `partition_id`.
#' @param call Caller environment for error reporting.
#' @return A tibble with `tract_id`, `zone_id`, `partition_id` (character),
#'   defaulting `partition_id` to `"all"` when absent from the input.
#' @keywords internal
.validate_zone_assignments <- function(zone_assignments, call = rlang::caller_env()) {
  if (!is.data.frame(zone_assignments)) {
    cli::cli_abort("{.arg zone_assignments} must be a data.frame.", call = call)
  }

  required <- c("tract_id", "zone_id")
  missing_cols <- setdiff(required, names(zone_assignments))
  if (length(missing_cols) > 0) {
    cli::cli_abort(
      "{.arg zone_assignments} is missing column{?s}: {.val {missing_cols}}.",
      call = call
    )
  }

  out <- tibble::tibble(
    tract_id = as.character(zone_assignments$tract_id),
    zone_id = as.character(zone_assignments$zone_id),
    partition_id = if ("partition_id" %in% names(zone_assignments)) {
      as.character(zone_assignments$partition_id)
    } else {
      "all"
    }
  )

  if (anyDuplicated(out$tract_id)) {
    dupes <- unique(out$tract_id[duplicated(out$tract_id)])
    cli::cli_abort(
      "{.field tract_id} values must be unique across zones: {.val {dupes}}.",
      call = call
    )
  }

  multi_partition <- out |>
    dplyr::distinct(.data$zone_id, .data$partition_id) |>
    dplyr::count(.data$zone_id) |>
    dplyr::filter(.data$n > 1) |>
    dplyr::pull(.data$zone_id)

  if (length(multi_partition) > 0) {
    cli::cli_abort(
      "Zone{?s} span more than one partition_id: {.val {multi_partition}}.",
      call = call
    )
  }

  out
}


#' Compute Symmetrized Total Distances Within a Zone
#'
#' @param tract_ids Character vector of tract IDs in one zone (length >= 2).
#' @param sparse_distances Tract-level sparse distance tibble.
#'
#' @return A named numeric vector (sum of symmetrized distances per
#'   candidate tract), or `NULL` if any pair has no distance in either
#'   direction.
#' @keywords internal
.zone_medoid_sums <- function(tract_ids, sparse_distances) {
  n <- length(tract_ids)

  within <- sparse_distances |>
    dplyr::filter(
      .data$origin_id %in% tract_ids,
      .data$destination_id %in% tract_ids,
      .data$origin_id != .data$destination_id
    ) |>
    dplyr::mutate(
      a = pmin(.data$origin_id, .data$destination_id),
      b = pmax(.data$origin_id, .data$destination_id)
    ) |>
    dplyr::summarise(
      distance = mean(.data$distance),
      .by = c("a", "b")
    )

  expected_pairs <- (n * (n - 1L)) %/% 2L
  if (nrow(within) < expected_pairs) {
    return(NULL)
  }

  sums <- stats::setNames(numeric(n), tract_ids)
  for (i in seq_len(nrow(within))) {
    sums[within$a[i]] <- sums[within$a[i]] + within$distance[i]
    sums[within$b[i]] <- sums[within$b[i]] + within$distance[i]
  }
  sums
}


#' Find Zone Centers by Minimum Total Distance (Medoid)
#'
#' For each zone, finds the tract minimizing the sum of (symmetrized)
#' distances to every other tract in the same zone. Zones with a single
#' tract trivially center on themselves.
#'
#' @param zone_assignments A validated tibble as returned by
#'   [.validate_zone_assignments()] (`tract_id`, `zone_id`, `partition_id`).
#' @param sparse_distances Tract-level sparse distance tibble (`origin_id`,
#'   `destination_id`, `distance`).
#' @param call Caller environment for error reporting.
#'
#' @return A tibble with `zone_id`, `partition_id`, `center_tract_id`,
#'   `n_tracts`.
#' @keywords internal
.find_zone_centers <- function(zone_assignments, sparse_distances,
                               call = rlang::caller_env()) {
  nested <- zone_assignments |>
    tidyr::nest(data = -c(zone_id, partition_id)) |>
    dplyr::mutate(
      n_tracts = purrr::map_int(.data$data, nrow),
      center_tract_id = purrr::map2_chr(
        .data$data, .data$n_tracts,
        \(d, n) {
          if (n == 1L) return(d$tract_id[1])
          sums <- .zone_medoid_sums(d$tract_id, sparse_distances)
          if (is.null(sums)) return(NA_character_)
          best <- min(sums)
          sort(names(sums)[sums == best])[1]
        }
      )
    )

  incomplete <- nested |>
    dplyr::filter(is.na(.data$center_tract_id)) |>
    dplyr::pull(.data$zone_id)

  if (length(incomplete) > 0) {
    cli::cli_abort(
      c(
        "Cannot find a center for zone{?s} with incomplete tract-pair distances: {.val {incomplete}}.",
        "i" = "Use {.code complete_distances = TRUE} with {.arg access_points} to fill missing pairs."
      ),
      call = call
    )
  }

  nested |>
    dplyr::select("zone_id", "partition_id", "center_tract_id", "n_tracts")
}

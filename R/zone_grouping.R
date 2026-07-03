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


#' Build a Sparse Zone-to-Zone Distance Table
#'
#' Wraps [.build_zone_dist_matrix()] (min pairwise tract distance) and
#' melts it into the sparse `origin_id`/`destination_id`/`distance` format
#' used elsewhere in the package, dropping self-pairs and any unreachable
#' (`Inf`) pairs.
#'
#' @param zone_ids Character vector of zone IDs.
#' @param zone_assignments Validated zone assignments (`tract_id`,
#'   `zone_id`).
#' @param sparse_distances Tract-level sparse distance tibble.
#'
#' @return A tibble with `origin_id`, `destination_id`, `distance`
#'   (zone IDs).
#' @keywords internal
.zone_dist_sparse <- function(zone_ids, zone_assignments, sparse_distances) {
  dist_mat <- .build_zone_dist_matrix(zone_ids, zone_assignments, sparse_distances)
  .matrix_to_sparse_dt(
    dist_mat, rownames(dist_mat), colnames(dist_mat), D_max = Inf
  ) |>
    dplyr::filter(is.finite(.data$distance))
}


#' Group Zones Within D_max Using the Reused Tract-Level p-Median Solver
#'
#' Treats zones as the atomic units of a second p-median solve: zone
#' centers become candidate "centers", and the same MILP machinery used
#' for tract-to-zone assignment now bounds every zone to within `D_max`
#' of its group's center zone.
#'
#' @param zone_ids Character vector of zone IDs in a single partition.
#' @param zone_assignments Validated zone assignments (`tract_id`,
#'   `zone_id`, `partition_id`) — used to compute the zone-to-zone
#'   distance matrix via minimum pairwise tract distance.
#' @param sparse_distances Tract-level sparse distance tibble.
#' @param D_max Numeric scalar. Grouping distance bound.
#' @param solver,max_time,rel_tol,verbose,use_cache Passed through to the
#'   reused solver, same meaning as in [surveyzones_build_zones()].
#'
#' @return A list with `assignments` (`zone_id`, `group_id`), `groups`
#'   (`group_id`, `center_zone_id`, `n_zones`, `diameter`), and
#'   `diagnostics`.
#' @keywords internal
.group_zones_pmedian <- function(zone_ids, zone_assignments, sparse_distances,
                                 D_max, solver, max_time, rel_tol, verbose,
                                 use_cache) {
  if (length(zone_ids) == 1L) {
    return(list(
      assignments = tibble::tibble(zone_id = zone_ids, group_id = zone_ids),
      groups = tibble::tibble(
        group_id = zone_ids, center_zone_id = zone_ids,
        n_zones = 1L, diameter = 0
      ),
      diagnostics = list(
        solver_status = "optimal", objective_value = 0,
        n_variables = 0L, solve_time = 0, optimality_gap = 0
      )
    ))
  }

  zone_sparse <- .zone_dist_sparse(zone_ids, zone_assignments, sparse_distances)

  synthetic_tracts <- tibble::tibble(
    tract_id = zone_ids,
    expected_service_time = 1
  )

  solve_single <- if (use_cache) {
    surveyzones_build_zones_single_mem
  } else {
    surveyzones_build_zones_single
  }

  result <- solve_single(
    full_sparse_distances = zone_sparse,
    tracts = synthetic_tracts,
    D_max = D_max,
    max_workload_per_zone = Inf,
    solver = solver,
    max_time = max_time,
    rel_tol = rel_tol,
    verbose = verbose,
    strategy = "auto"
  )

  list(
    assignments = tibble::tibble(
      zone_id = result$assignments$tract_id,
      group_id = result$assignments$zone_id
    ),
    groups = tibble::tibble(
      group_id = result$zones$zone_id,
      center_zone_id = result$zones$center_tract_id,
      n_zones = result$zones$n_tracts,
      diameter = result$zones$diameter
    ),
    diagnostics = result$diagnostics
  )
}

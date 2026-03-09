#' Sequence All Zones in a Plan
#'
#' Computes a visit order for every zone's tracts.
#'
#' @param plan A `surveyzones_plan` object.
#' @param sparse_distances Sparse distance table (as returned by
#'   [surveyzones_compute_sparse_distances()]).
#' @param method Character scalar passed to [seriation::seriate()]
#'   (e.g., `"TSP"`, `"SPIN_NH"`, `"Spectral"`, `"OLO"`).
#'   Default `"TSP"`.
#' @param control Named list of method-specific parameters passed to
#'   [seriation::seriate()].  For example,
#'   `control = list(sigma = c(10, 7, 5, 3), step = 10)` for `"SPIN_NH"`.
#'   Default `NULL`.
#'
#' @return The same `plan` with `plan$sequence` populated — a tibble
#'   with columns `zone_id`, `tract_id`, `visit_order`.
#'
#' @export
surveyzones_sequence <- function(
  plan,
  sparse_distances,
  method = "TSP",
  control = NULL
) {
  if (!inherits(plan, "surveyzones_plan")) {
    cli::cli_abort("{.arg plan} must be a {.cls surveyzones_plan} object.")
  }

  plan$sequence <- plan$assignments |>
    tidyr::nest(data = -zone_id) |>
    dplyr::mutate(
      ordered = purrr::map(
        data,
        \(z) {
          surveyzones_sequence_tracts(
            tract_ids = z$tract_id,
            sparse_distances = sparse_distances,
            method = method,
            control = control
          )
        }
      )
    ) |>
    dplyr::select(-data) |>
    tidyr::unnest_longer(ordered, values_to = "tract_id") |>
    dplyr::group_by(zone_id) |>
    dplyr::mutate(visit_order = dplyr::row_number()) |>
    dplyr::ungroup() |>
    dplyr::select(zone_id, tract_id, visit_order)

  plan
}


#' Seriate a Set of IDs into 1D Order
#'
#' Builds a symmetric distance matrix from a sparse distance table and
#' applies [seriation::seriate()] to find a 1D ordering.
#'
#' @param ids Character vector of IDs to order.
#' @param sparse_distances Sparse distance tibble.
#' @param method Character scalar passed to [seriation::seriate()].
#' @param control Named list of method-specific parameters passed to
#'   [seriation::seriate()].
#'
#' @return Character vector of IDs in seriation order.
#' @keywords internal
.seriate_1d <- function(ids, sparse_distances, method = "TSP", control = NULL) {
  n <- length(ids)
  if (n <= 1) {
    return(ids)
  }

  # Build complete symmetric distance matrix
  mat <- matrix(0, nrow = n, ncol = n, dimnames = list(ids, ids))

  dt <- sparse_distances |>
    dplyr::filter(
      origin_id %in% ids,
      destination_id %in% ids,
      origin_id != destination_id
    )

  # Check completeness: need n*(n-1) directed pairs
  expected_pairs <- n * (n - 1)
  unique_pairs <- nrow(dplyr::distinct(dt, origin_id, destination_id))
  if (unique_pairs < expected_pairs) {
    n_missing <- expected_pairs - unique_pairs
    cli::cli_abort(c(
      "Distance matrix is incomplete ({n_missing} missing pair{?s}).",
      "i" = "Use {.fun surveyzones_complete_distances} to fill missing pairs before sequencing."
    ))
  }

  # Fill matrix, symmetrise by averaging d(a,b) and d(b,a)
  if (nrow(dt) > 0) {
    mat[cbind(
      as.character(dt$origin_id),
      as.character(dt$destination_id)
    )] <- dt$distance
  }
  mat <- (mat + t(mat)) / 2

  # Apply SPIN_NH defaults:
  if (method == "SPIN_NH" && is.null(control)) {
    control <- list(sigma = 1)
  }

  # Seriate
  o <- seriation::seriate(
    stats::as.dist(mat),
    method = method,
    control = control
  )
  perm <- seriation::get_order(o)
  ids[perm]
}


#' Build a Zone-to-Zone Distance Matrix Using Min Pairwise Distances
#'
#' Computes `d(zone_A, zone_B) = min(d(tract_i, tract_j))` for all
#' `i in zone_A, j in zone_B` using the sparse distance table.
#'
#' @param zone_ids Character vector of zone IDs.
#' @param assignments Tibble with `tract_id` and `zone_id` columns.
#' @param sparse_distances Sparse distance tibble.
#'
#' @return A symmetric numeric matrix with zone_ids as row/col names.
#'   Diagonal is 0. Missing zone pairs get `Inf`.
#' @keywords internal
.build_zone_dist_matrix <- function(zone_ids, assignments, sparse_distances) {
  n <- length(zone_ids)
  mat <- matrix(Inf, nrow = n, ncol = n, dimnames = list(zone_ids, zone_ids))
  diag(mat) <- 0

  if (n <= 1) return(mat)

  # Map tract_id -> zone_id
  tract_zone <- assignments |>
    dplyr::select("tract_id", "zone_id")

  # Join zone_ids onto sparse distances
  cross <- sparse_distances |>
    dplyr::inner_join(tract_zone, by = dplyr::join_by(origin_id == tract_id)) |>
    dplyr::rename(zone_origin = "zone_id") |>
    dplyr::inner_join(tract_zone, by = dplyr::join_by(destination_id == tract_id)) |>
    dplyr::rename(zone_dest = "zone_id") |>
    dplyr::filter(
      zone_origin != zone_dest,
      zone_origin %in% zone_ids,
      zone_dest %in% zone_ids
    )

  if (nrow(cross) == 0) return(mat)

  # Min distance per zone pair
  zone_dists <- cross |>
    dplyr::summarise(
      min_dist = min(distance),
      .by = c(zone_origin, zone_dest)
    )

  # Fill matrix
  mat[cbind(zone_dists$zone_origin, zone_dists$zone_dest)] <- zone_dists$min_dist

  # Symmetrise: take min of d(A,B) and d(B,A)
  mat <- pmin(mat, t(mat))

  mat
}

#' Sequence Zones Within Each Partition
#'
#' Given a solved plan, finds a spatial ordering for the zones within each
#' partition using TSP seriation on minimum pairwise tract distances.
#' Groups are assigned by scanning consecutive distances in the TSP path
#' and splitting wherever the distance exceeds `threshold`.
#'
#' @param plan A `surveyzones_plan` object.
#' @param sparse_distances Sparse distance table (as returned by
#'   [surveyzones_compute_sparse_distances()]).
#' @param access_points An sf object with POINT geometries and a
#'   `tract_id` column.  Used to fill missing pairwise distances via
#'   haversine ([surveyzones_complete_distances()]) so that every
#'   cross-zone tract pair has a distance.
#'   Defaults to `plan$access_points` when `NULL`.  Required when
#'   `complete_distances = TRUE`.
#' @param speed_kmh Numeric scalar.  Assumed travel speed for haversine
#'   fill-in (km/h).  Default `0.1` — intentionally harsh because missing
#'   pairs likely represent real barriers (rivers, mountains).
#' @param complete_distances Logical scalar.  When `TRUE` (default),
#'   missing pairwise distances are filled using haversine via
#'   [surveyzones_complete_distances()].  Set to `FALSE` to skip
#'   imputation (e.g., when distances are already complete).
#' @param threshold Numeric scalar.  Maximum consecutive distance (in the
#'   same units as `sparse_distances`, typically minutes) between two zones
#'   in the TSP path before starting a new group.  Default `10`.
#' @param by_partition Logical scalar.  When `TRUE` (default), zones are
#'   sequenced independently within each partition.  When `FALSE`, all zones
#'   are sequenced together ignoring partition boundaries.
#'
#' @return The same `plan` with `plan$zone_sequence` populated — a tibble
#'   with columns `partition_id`, `zone_id`, `zone_order`, `group_id`.
#'
#' @export
surveyzones_sequence_zones <- function(
  plan,
  sparse_distances,
  access_points = NULL,
  speed_kmh = 0.1,
  complete_distances = TRUE,
  threshold = 10,
  by_partition = TRUE
) {
  if (!inherits(plan, "surveyzones_plan")) {
    cli::cli_abort("{.arg plan} must be a {.cls surveyzones_plan} object.")
  }

  access_points <- access_points %||% plan$access_points

  # Complete missing distances using haversine so every cross-zone pair exists
  if (complete_distances) {
    if (is.null(access_points)) {
      cli::cli_abort(c(
        "{.arg access_points} is required when {.arg complete_distances} is {.val TRUE}.",
        "i" = "Either pass it explicitly, store it in the plan via {.fun surveyzones_build_zones}, or set {.code complete_distances = FALSE}."
      ))
    }
    sparse_distances <- surveyzones_complete_distances(
      sparse_distances, access_points, speed_kmh = speed_kmh
    )
  }

  # Determine grouping column
  group_col <- if (by_partition) "partition_id" else ".group"
  zones <- plan$zones
  if (!by_partition) {
    zones <- dplyr::mutate(zones, .group = "all")
  }

  zone_lookup <- plan$zones |>
    dplyr::select("zone_id", "partition_id")

  plan$zone_sequence <- zones |>
    tidyr::nest(data = -dplyr::all_of(group_col)) |>
    dplyr::mutate(
      ordered = purrr::map(data, \(z) {
        zone_ids <- z$zone_id
        n <- length(zone_ids)
        if (n <= 1) {
          return(tibble::tibble(zone_id = zone_ids, group_id = 1L))
        }

        # Build min-pairwise zone distance matrix
        mat <- .build_zone_dist_matrix(
          zone_ids, plan$assignments, sparse_distances
        )

        # TSP-seriate all zones into a single optimal path
        zone_sparse <- .matrix_to_sparse(mat)
        ordered_ids <- .seriate_1d(zone_ids, zone_sparse, method = "TSP")

        # Scan consecutive distances and split at gaps > threshold
        consecutive_dists <- vapply(
          seq_len(n - 1),
          \(i) mat[ordered_ids[i], ordered_ids[i + 1]],
          numeric(1)
        )
        group_id <- c(1L, cumsum(consecutive_dists > threshold) + 1L)

        tibble::tibble(zone_id = ordered_ids, group_id = group_id)
      })
    ) |>
    dplyr::select(-data) |>
    tidyr::unnest(ordered) |>
    dplyr::group_by(dplyr::across(dplyr::all_of(group_col))) |>
    dplyr::mutate(zone_order = dplyr::row_number()) |>
    dplyr::ungroup()

  # When by_partition = FALSE, restore partition_id from the zone lookup
  if (!by_partition) {
    plan$zone_sequence <- plan$zone_sequence |>
      dplyr::left_join(zone_lookup, by = "zone_id") |>
      dplyr::select("partition_id", "zone_id", "zone_order", "group_id")
  } else {
    plan$zone_sequence <- plan$zone_sequence |>
      dplyr::select("partition_id", "zone_id", "zone_order", "group_id")
  }

  # Orient tract sequences so each zone's first tract is closest to the

  # previous zone's last tract
  if (!is.null(plan$sequence)) {
    plan$sequence <- .orient_tract_sequences(
      plan$sequence, plan$zone_sequence, sparse_distances
    )
  }

  plan
}


#' Orient Tract Sequences Based on Zone Visit Order
#'
#' For each zone after the first, checks whether the tract sequence should
#' be reversed so that the entry tract (first in visit order) is the one
#' closest to the exit tract (last in visit order) of the previous zone.
#'
#' @param sequence Tibble with `zone_id`, `tract_id`, `visit_order`.
#' @param zone_sequence Tibble with `partition_id`, `zone_id`, `zone_order`.
#' @param sparse_distances Sparse distance tibble.
#'
#' @return Updated `sequence` tibble with potentially reversed tract orders.
#' @keywords internal
.orient_tract_sequences <- function(sequence, zone_sequence, sparse_distances) {
  # Process each partition independently
  partitions <- split(zone_sequence, zone_sequence$partition_id)

  for (part in partitions) {
    ordered_zones <- part$zone_id[order(part$zone_order)]
    if (length(ordered_zones) <= 1) next

    for (i in seq(2, length(ordered_zones))) {
      prev_zone <- ordered_zones[i - 1]
      curr_zone <- ordered_zones[i]

      # Get the last tract of the previous zone (exit point)
      prev_tracts <- sequence |>
        dplyr::filter(zone_id == prev_zone) |>
        dplyr::arrange(visit_order)
      exit_tract <- prev_tracts$tract_id[nrow(prev_tracts)]

      # Get the first and last tract of the current zone
      curr_tracts <- sequence |>
        dplyr::filter(zone_id == curr_zone) |>
        dplyr::arrange(visit_order)
      if (nrow(curr_tracts) <= 1) next

      first_tract <- curr_tracts$tract_id[1]
      last_tract <- curr_tracts$tract_id[nrow(curr_tracts)]

      # Look up distances from exit to first vs exit to last
      d_to_first <- sparse_distances |>
        dplyr::filter(origin_id == exit_tract, destination_id == first_tract) |>
        dplyr::pull(distance)
      d_to_last <- sparse_distances |>
        dplyr::filter(origin_id == exit_tract, destination_id == last_tract) |>
        dplyr::pull(distance)

      # If closer to last tract, reverse the sequence
      if (length(d_to_first) > 0 && length(d_to_last) > 0 &&
          d_to_last[1] < d_to_first[1]) {
        n <- nrow(curr_tracts)
        sequence$visit_order[sequence$zone_id == curr_zone] <-
          rev(sequence$visit_order[sequence$zone_id == curr_zone])
      }
    }
  }

  sequence
}


#' Convert a Symmetric Distance Matrix to Sparse Format
#'
#' @param mat Symmetric numeric matrix with IDs as row/col names.
#' @return A tibble with columns `origin_id`, `destination_id`, `distance`.
#' @keywords internal
.matrix_to_sparse <- function(mat) {
  ids <- rownames(mat)
  n <- length(ids)
  if (n <= 1) {
    return(tibble::tibble(
      origin_id = character(0),
      destination_id = character(0),
      distance = numeric(0)
    ))
  }
  pairs <- tidyr::expand_grid(origin_id = ids, destination_id = ids) |>
    dplyr::filter(origin_id != destination_id) |>
    dplyr::mutate(
      distance = mat[cbind(origin_id, destination_id)]
    )
  pairs
}

#' Sequence Tracts Within a Single Zone
#'
#' Orders tracts spatially using [seriation::seriate()] for 1D ordering.
#'
#' @param tract_ids Character vector of tract identifiers in this zone.
#' @param sparse_distances Sparse distance table.
#' @param method Character scalar passed to [seriation::seriate()]
#'   (e.g., `"TSP"`, `"SPIN_NH"`, `"Spectral"`, `"OLO"`).
#'   Default `"TSP"`.
#' @param control Named list of method-specific parameters passed to
#'   [seriation::seriate()].  Default `NULL`.
#'
#' @return Character vector of tract_ids in visit order.
#'
#' @export
surveyzones_sequence_tracts <- function(
  tract_ids,
  sparse_distances,
  method = "TSP",
  control = NULL
) {
  n <- length(tract_ids)
  if (n <= 2) {
    return(tract_ids)
  }

  .seriate_1d(tract_ids, sparse_distances, method = method, control = control)
}

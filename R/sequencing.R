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

#' Flood-Fill Zones into Threshold-Connected Groups
#'
#' Partitions zones into groups where every zone in a group is reachable
#' from every other zone in the same group via a chain of pairwise
#' distances <= `threshold`.
#'
#' @param mat Symmetric numeric distance matrix with zone IDs as row/col names.
#' @param threshold Numeric scalar.  Maximum distance for two zones to be
#'   considered connected.
#'
#' @return A list of character vectors, each containing the zone IDs in one
#'   group.
#' @keywords internal
.flood_fill_groups <- function(mat, threshold) {
  zone_ids <- rownames(mat)
  n <- length(zone_ids)
  visited <- rep(FALSE, n)
  names(visited) <- zone_ids
  groups <- list()

  for (seed in zone_ids) {
    if (visited[seed]) next

    # BFS from seed
    queue <- seed
    group <- character(0)
    while (length(queue) > 0) {
      current <- queue[1]
      queue <- queue[-1]
      if (visited[current]) next
      visited[current] <- TRUE
      group <- c(group, current)

      # Find unvisited neighbors within threshold
      neighbors <- zone_ids[!visited & mat[current, ] <= threshold & mat[current, ] > 0]
      queue <- c(queue, neighbors)
    }

    groups <- c(groups, list(group))
  }

  groups
}

#' Sequence Zones Within Each Partition
#'
#' Given a solved plan, finds a spatial ordering for the zones within each
#' partition using a flood-fill + TSP approach.  Zones within `threshold`
#' distance of each other are grouped together, then groups and zones within
#' groups are ordered via TSP seriation on minimum pairwise tract distances.
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
#' @param threshold Numeric scalar.  Maximum distance (in the same units as
#'   `sparse_distances`, typically minutes) for two zones to be grouped
#'   together.  Default `10`.
#' @param by_partition Logical scalar.  When `TRUE` (default), zones are
#'   sequenced independently within each partition.  When `FALSE`, all zones
#'   are sequenced together ignoring partition boundaries.
#'
#' @return The same `plan` with `plan$zone_sequence` populated — a tibble
#'   with columns `partition_id`, `zone_id`, `zone_order`.
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
        if (n <= 1) return(zone_ids)

        # Build min-pairwise zone distance matrix
        mat <- .build_zone_dist_matrix(
          zone_ids, plan$assignments, sparse_distances
        )

        # Flood-fill into threshold-connected groups
        groups <- .flood_fill_groups(mat, threshold)

        # Convert zone distance matrix to sparse format for .seriate_1d()
        zone_sparse <- .matrix_to_sparse(mat)

        if (length(groups) == 1) {
          # Single group: just TSP-seriate all zones
          return(.seriate_1d(zone_ids, zone_sparse, method = "TSP"))
        }

        # Seriate within each group
        ordered_groups <- purrr::map(groups, \(g) {
          if (length(g) <= 2) return(g)
          .seriate_1d(g, zone_sparse, method = "TSP")
        })

        # Build group-to-group distance matrix
        ng <- length(ordered_groups)
        group_names <- paste0("G", seq_len(ng))
        group_mat <- matrix(Inf, ng, ng, dimnames = list(group_names, group_names))
        diag(group_mat) <- 0

        for (i in seq_len(ng - 1)) {
          for (j in (i + 1):ng) {
            d <- min(mat[groups[[i]], groups[[j]], drop = FALSE])
            group_mat[i, j] <- d
            group_mat[j, i] <- d
          }
        }

        # Seriate groups
        if (ng <= 2) {
          group_order <- seq_len(ng)
        } else {
          group_dist <- stats::as.dist(group_mat)
          o <- seriation::seriate(group_dist, method = "TSP")
          group_order <- seriation::get_order(o)
        }

        # Concatenate: group order × within-group order
        unlist(ordered_groups[group_order])
      })
    ) |>
    dplyr::select(-data) |>
    tidyr::unnest_longer(ordered, values_to = "zone_id") |>
    dplyr::group_by(dplyr::across(dplyr::all_of(group_col))) |>
    dplyr::mutate(zone_order = dplyr::row_number()) |>
    dplyr::ungroup()

  # When by_partition = FALSE, restore partition_id from the zone lookup
  if (!by_partition) {
    plan$zone_sequence <- plan$zone_sequence |>
      dplyr::left_join(zone_lookup, by = "zone_id") |>
      dplyr::select("partition_id", "zone_id", "zone_order")
  } else {
    plan$zone_sequence <- plan$zone_sequence |>
      dplyr::select("partition_id", "zone_id", "zone_order")
  }

  plan
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

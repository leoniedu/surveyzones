#' Sequence Zones and Tracts in a Plan
#'
#' Computes a spatial visit order at two levels: (1) zones within each
#' partition via TSP on minimum pairwise tract distances, and (2) tracts
#' within each zone via seriation.  Tract sequences are oriented so that
#' each zone's entry tract is closest to the previous zone's exit tract.
#'
#' @param plan A `surveyzones_plan` object.
#' @param sparse_distances Sparse distance table (as returned by
#'   [surveyzones_compute_sparse_distances()]).
#' @param method Character scalar passed to [seriation::seriate()]
#'   for tract-level ordering
#'   (e.g., `"TSP"`, `"SPIN_NH"`, `"Spectral"`, `"OLO"`).
#'   Default `"TSP"`.
#' @param control Named list of method-specific parameters passed to
#'   [seriation::seriate()] (and through to [TSP::solve_TSP()] when
#'   `method = "TSP"`).  Default `NULL` (seriation's default heuristic).
#'   Example: `list(method = "nn", rep = 20)` for nearest-neighbour.
#'   When `method` is not `"TSP"`, TSP-specific control keys are
#'   silently stripped.
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
#' @param by_partition Logical scalar.  When `TRUE` (default), zones are
#'   sequenced independently within each partition.  When `FALSE`, all zones
#'   are sequenced together ignoring partition boundaries.
#'
#' @return The same `plan` with two new components:
#'   - `plan$zone_sequence`: tibble with `partition_id`, `zone_id`,
#'     `zone_order`, `group_id`.
#'   - `plan$sequence`: tibble with `zone_id`, `tract_id`, `visit_order`.
#'
#' @export
surveyzones_sequence <- function(
  plan,
  sparse_distances,
  method = "TSP",
  control = NULL,
  access_points = NULL,
  speed_kmh = 0.1,
  complete_distances = TRUE,
  by_partition = TRUE
) {
  if (!inherits(plan, "surveyzones_plan")) {
    cli::cli_abort("{.arg plan} must be a {.cls surveyzones_plan} object.")
  }

  access_points <- access_points %||% plan$access_points

  # Complete missing distances using haversine
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

  plan$zone_sequence <- .sequence_zones(
    plan, sparse_distances, control, by_partition
  )

  plan$sequence <- .sequence_tracts_with_orientation(
    plan, sparse_distances, method, control
  )

  .rename_zones(plan)
}


#' Sequence Zones Within Partitions
#'
#' TSP-seriates zones using min-pairwise tract distances, then splits
#' using the solver-defined groups.
#'
#' @param plan A `surveyzones_plan` object.
#' @param sparse_distances Sparse distance tibble.
#' @param control Control list passed to [seriation::seriate()].
#' @param by_partition Whether to sequence within partitions.
#'
#' @return Tibble with `partition_id`, `zone_id`, `zone_order`, `group_id`.
#' @keywords internal
.sequence_zones <- function(plan, sparse_distances, control, by_partition) {
  zones <- plan$zones

  # Save original partition_id for output
  zone_lookup <- zones |> dplyr::select("zone_id", orig_partition = "partition_id")

  # When by_partition = FALSE, treat all zones as one partition
  if (!by_partition) {
    zones$partition_id <- "all"
  }

  # Step 1: Order zones within each (partition_id, group_id) using TSP
  within_group <- zones |>
    tidyr::nest(data = -dplyr::all_of(c("partition_id", "group_id"))) |>
    dplyr::mutate(
      ordered = purrr::map(data, \(z) {
        zone_ids <- z$zone_id
        n <- length(zone_ids)
        if (n <= 1) {
          return(tibble::tibble(zone_id = zone_ids, within_order = 1L))
        }

        mat <- .build_zone_dist_matrix(
          zone_ids, plan$assignments, sparse_distances
        )

        o <- .seriate_quietly(stats::as.dist(mat), method = "TSP",
                              control = control)
        ordered_ids <- zone_ids[seriation::get_order(o)]

        tibble::tibble(
          zone_id = ordered_ids,
          within_order = seq_along(ordered_ids)
        )
      })
    ) |>
    dplyr::select(-data) |>
    tidyr::unnest(ordered)

  # Step 2: Order groups within each partition using TSP on group medoids
  # (pick first zone in each group as representative)
  group_order <- within_group |>
    dplyr::filter(.data$within_order == 1L) |>
    tidyr::nest(gdata = -"partition_id") |>
    dplyr::mutate(
      gordered = purrr::map(gdata, \(g) {
        group_ids <- g$group_id
        rep_zone_ids <- g$zone_id
        n <- length(group_ids)
        if (n <= 1) {
          return(tibble::tibble(group_id = group_ids, group_order = 1L))
        }
        mat <- .build_zone_dist_matrix(
          rep_zone_ids, plan$assignments, sparse_distances
        )
        o <- .seriate_quietly(stats::as.dist(mat), method = "TSP",
                              control = control)
        tibble::tibble(
          group_id = group_ids[seriation::get_order(o)],
          group_order = seq_len(n)
        )
      })
    ) |>
    dplyr::select(-gdata) |>
    tidyr::unnest(gordered)

  # Step 3: Combine and compute global zone_order within partition
  result <- within_group |>
    dplyr::left_join(group_order, by = c("partition_id", "group_id")) |>
    dplyr::arrange(.data$partition_id, .data$group_order, .data$within_order) |>
    dplyr::mutate(
      zone_order = dplyr::row_number(),
      .by = "partition_id"
    ) |>
    dplyr::select("partition_id", "zone_id", "zone_order", "group_id")

  # Restore original partition_id when by_partition = FALSE
  if (!by_partition) {
    result <- result |>
      dplyr::select(-"partition_id") |>
      dplyr::left_join(zone_lookup, by = "zone_id") |>
      dplyr::rename(partition_id = "orig_partition") |>
      dplyr::select("partition_id", "zone_id", "zone_order", "group_id")
  }

  result
}


#' Sequence and Orient Tracts Within Each Zone
#'
#' Walks zones in visit order, seriates tracts within each zone, and
#' orients each zone's tract sequence so the entry tract is closest
#' to the previous zone's exit tract.
#'
#' @param plan A `surveyzones_plan` with `zone_sequence` already set.
#' @param sparse_distances Sparse distance tibble.
#' @param method Seriation method for tracts.
#' @param control Control list passed to [seriation::seriate()].
#'
#' @return Tibble with `zone_id`, `tract_id`, `visit_order`.
#' @keywords internal
.sequence_tracts_with_orientation <- function(plan, sparse_distances,
                                              method, control) {
  # Build O(1) distance lookup
  dist_key <- paste(
    sparse_distances$origin_id, sparse_distances$destination_id,
    sep = "|"
  )
  dist_val <- sparse_distances$distance
  names(dist_val) <- dist_key

  # Group tracts by zone
  tracts_by_zone <- split(
    plan$assignments$tract_id,
    plan$assignments$zone_id
  )

  # Walk zones in visit order per partition; orient as we go
  zs <- plan$zone_sequence |> dplyr::arrange(partition_id, zone_order)
  sequences <- vector("list", nrow(zs))
  prev_exit <- character(0)  # named by partition_id

  for (j in seq_len(nrow(zs))) {
    zid <- zs$zone_id[j]
    pid <- zs$partition_id[j]
    tract_ids <- tracts_by_zone[[zid]]

    ordered_tracts <- surveyzones_sequence_tracts(
      tract_ids, sparse_distances, method = method, control = control
    )

    # Orient: if previous zone exists in this partition and ≥2 tracts,
    # check whether reversing brings entry closer to previous exit
    exit <- prev_exit[pid]
    if (!is.na(exit) && length(ordered_tracts) > 1) {
      first <- ordered_tracts[1]
      last <- ordered_tracts[length(ordered_tracts)]
      d_first <- dist_val[paste(exit, first, sep = "|")]
      d_last <- dist_val[paste(exit, last, sep = "|")]
      if (!is.na(d_first) && !is.na(d_last) && d_last < d_first) {
        ordered_tracts <- rev(ordered_tracts)
      }
    }

    prev_exit[pid] <- ordered_tracts[length(ordered_tracts)]

    sequences[[j]] <- tibble::tibble(
      zone_id = zid,
      tract_id = ordered_tracts,
      visit_order = seq_along(ordered_tracts)
    )
  }

  dplyr::bind_rows(sequences)
}


#' Call seriation::seriate() Suppressing Only TSP Control Warnings
#'
#' @param d A `dist` object.
#' @param method Seriation method.
#' @param control Control list.
#'
#' @return A seriation order object.
#' @keywords internal
.seriate_quietly <- function(d, method, control) {
  withCallingHandlers(
    seriation::seriate(d, method = method, control = control),
    warning = function(w) {
      if (grepl("control parameter|Unknown parameter|sequentially", conditionMessage(w))) {
        invokeRestart("muffleWarning")
      }
    }
  )
}


#' Rename zone IDs to sequential {partition_id}_{group_id}.{pos} format
#'
#' Applied at the end of [surveyzones_sequence()] after zone_sequence is
#' computed.  Builds an old->new lookup and applies it to all four plan tables.
#' Also drops `group_id` from `assignments` (character phase-1 center tract ID,
#' redundant with `zones$center_tract_id`). Note: `zone_sequence$group_id` (integer
#' travel-group counter) is used in the new ID format — these are different columns.
#'
#' @param plan A `surveyzones_plan` with `zone_sequence` already set.
#' @return The plan with renamed zone_ids in all four tables.
#' @keywords internal
.rename_zones <- function(plan) {
  # Build lookup: assign sequential group index within partition,
  # then pos = rank within (partition_id, group_idx) ordered by zone_order
  group_idx_lookup <- plan$zone_sequence |>
    dplyr::distinct(.data$partition_id, .data$group_id) |>
    dplyr::mutate(
      group_idx = dplyr::row_number(),
      .by = "partition_id"
    )

  lookup <- plan$zone_sequence |>
    dplyr::left_join(group_idx_lookup, by = c("partition_id", "group_id")) |>
    dplyr::arrange(.data$partition_id, .data$group_idx, .data$zone_order) |>
    dplyr::mutate(
      pos = rank(.data$zone_order, ties.method = "first"),
      .by = c("partition_id", "group_idx")
    ) |>
    dplyr::mutate(
      new_zone_id = paste0(
        .data$partition_id, "_",
        .data$group_idx, ".",
        sprintf("%03d", .data$pos)
      )
    ) |>
    dplyr::select("zone_id", "new_zone_id")

  id_map <- stats::setNames(lookup$new_zone_id, lookup$zone_id)

  missing <- setdiff(
    c(plan$zone_sequence$zone_id, plan$zones$zone_id,
      plan$assignments$zone_id, plan$sequence$zone_id),
    names(id_map)
  )
  if (length(missing) > 0) {
    cli::cli_abort("zone_ids not found in lookup: {.val {missing}}")
  }

  # Apply to all four tables
  plan$zone_sequence$zone_id <- unname(id_map[plan$zone_sequence$zone_id])
  plan$zones$zone_id         <- unname(id_map[plan$zones$zone_id])
  plan$assignments$zone_id   <- unname(id_map[plan$assignments$zone_id])
  plan$sequence$zone_id      <- unname(id_map[plan$sequence$zone_id])

  # Drop character group_id from assignments (redundant with zones$center_tract_id)
  plan$assignments <- dplyr::select(plan$assignments, -dplyr::any_of("group_id"))

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
#'   [seriation::seriate()].  Default `NULL`.
#'
#' @return Character vector of IDs in seriation order.
#' @keywords internal
.seriate_1d <- function(ids, sparse_distances, method = "TSP",
                        control = NULL) {
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

  # Strip TSP-specific keys when using non-TSP methods
  if (method != "TSP" && !is.null(control)) {
    tsp_keys <- c("method", "rep", "two_opt", "start")
    control <- control[setdiff(names(control), tsp_keys)]
    if (length(control) == 0) control <- NULL
  }

  # Apply SPIN_NH defaults:
  if (method == "SPIN_NH" && is.null(control)) {
    control <- list(sigma = 1)
  }

  o <- .seriate_quietly(d = stats::as.dist(mat), method = method,
                        control = control)
  ids[seriation::get_order(o)]
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

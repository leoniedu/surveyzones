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


#' Reorient a Zone/Group Visit Order Toward a Partition's Starting Point
#'
#' Compares the endpoint tracts of an ordered ID sequence (each ID mapped
#' to its representative tract) to a partition's starting-point
#' coordinate, and reverses the sequence when its *last* element is
#' closer to the starting point than its *first* element -- so the
#' sequence starts at the (endpoint-optimal) closer end. TSP paths have
#' exactly two valid global orientations; this picks the better one
#' without re-solving the TSP with a fixed start node.
#'
#' @param order_ids Character vector, already TSP-ordered (zone IDs or
#'   group center-zone IDs).
#' @param zone_centers Named character vector: representative `tract_id`
#'   for each ID in `order_ids` (names = the IDs).
#' @param starting_point A one-row sf POINT object (the partition's
#'   starting point), or `NULL` to skip anchoring.
#' @param access_points An sf object with POINT geometries and a
#'   `tract_id` column, covering every tract in `zone_centers`.
#'
#' @return `order_ids`, possibly reversed.
#' @keywords internal
.anchor_partition_order <- function(order_ids, zone_centers, starting_point,
                                    access_points) {
  if (length(order_ids) <= 1L || is.null(starting_point)) return(order_ids)

  first_tract <- zone_centers[[order_ids[1]]]
  last_tract <- zone_centers[[order_ids[length(order_ids)]]]

  first_pt <- access_points[access_points$tract_id == first_tract, ]
  last_pt <- access_points[access_points$tract_id == last_tract, ]

  d_first <- as.numeric(sf::st_distance(starting_point, first_pt))
  d_last <- as.numeric(sf::st_distance(starting_point, last_pt))

  if (d_last < d_first) rev(order_ids) else order_ids
}


#' Assign ggg.zzz Labels via Two Independent Partition-Wide TSPs
#'
#' Seriates group centers (for `ggg`) and all zones (for `zzz`)
#' independently within one partition -- the two orders do not need to
#' agree, so a group's members may be non-contiguous in the `zzz`
#' sequence. When a starting point is supplied, both orders are anchored
#' via [.anchor_partition_order()] so the first zone and first group are
#' each the one closest to it.
#'
#' @param zone_ids Character vector of all zone IDs in one partition.
#' @param groups Tibble for this partition with `group_id`,
#'   `center_zone_id` (as returned by [.group_zones_pmedian()]).
#' @param zone_assignments Validated zone assignments.
#' @param sparse_distances Tract-level sparse distance tibble.
#' @param zone_centers Named character vector (`zone_id` -> representative
#'   `tract_id`, scoped to this partition). Only used when
#'   `starting_point` is supplied.
#' @param starting_point A one-row sf POINT object (this partition's
#'   starting point), or `NULL` (default) to skip anchoring.
#' @param access_points An sf object with POINT geometries and a
#'   `tract_id` column. Required when `starting_point` is supplied.
#'
#' @return A list with `zone_labels` (`zone_id`, `zzz`) and
#'   `group_labels` (`group_id`, `ggg`).
#' @keywords internal
.label_zone_groups <- function(zone_ids, groups, zone_assignments,
                               sparse_distances, zone_centers = NULL,
                               starting_point = NULL, access_points = NULL) {
  zone_dists <- .zone_dist_sparse(zone_ids, zone_assignments, sparse_distances)

  zone_order <- .seriate_1d(zone_ids, zone_dists)
  group_order <- .seriate_1d(groups$center_zone_id, zone_dists)

  zone_order <- .anchor_partition_order(
    zone_order, zone_centers, starting_point, access_points
  )
  group_order <- .anchor_partition_order(
    group_order, zone_centers, starting_point, access_points
  )

  zone_labels <- tibble::tibble(
    zone_id = zone_order,
    zzz = sprintf("%03d", seq_along(zone_order))
  )

  group_labels <- tibble::tibble(
    center_zone_id = group_order,
    ggg = sprintf("%03d", seq_along(group_order))
  ) |>
    dplyr::left_join(
      groups |> dplyr::select("group_id", "center_zone_id"),
      by = "center_zone_id"
    ) |>
    dplyr::select("group_id", "ggg")

  list(zone_labels = zone_labels, group_labels = group_labels)
}


#' Sequence and Orient Tracts Along a Given Zone Visit Order
#'
#' Walks zones in the supplied order, seriates each zone's tracts, and
#' orients each zone's tract sequence so its entry tract is closest to the
#' previous zone's exit tract. Mirrors
#' `.sequence_tracts_with_orientation()` in `R/sequencing.R`, but driven
#' by an externally supplied zone order (here, the `zzz` TSP order from
#' [.label_zone_groups()]) instead of a group+within-group order.
#'
#' @param zone_order Character vector of zone IDs, in visit order
#'   (one partition).
#' @param zone_assignments Validated zone assignments (`tract_id`,
#'   `zone_id`).
#' @param sparse_distances Tract-level sparse distance tibble.
#' @param method,control Passed to [surveyzones_sequence_tracts()].
#'
#' @return A tibble with `zone_id`, `tract_id`, `visit_order`.
#' @keywords internal
.sequence_group_tracts <- function(zone_order, zone_assignments,
                                   sparse_distances, method = "TSP",
                                   control = NULL) {
  dist_key <- paste(
    sparse_distances$origin_id, sparse_distances$destination_id, sep = "|"
  )
  dist_val <- stats::setNames(sparse_distances$distance, dist_key)

  tracts_by_zone <- split(
    zone_assignments$tract_id, zone_assignments$zone_id
  )

  sequences <- vector("list", length(zone_order))
  prev_exit <- NA_character_

  for (j in seq_along(zone_order)) {
    zid <- zone_order[j]
    tract_ids <- tracts_by_zone[[zid]]

    ordered_tracts <- surveyzones_sequence_tracts(
      tract_ids, sparse_distances, method = method, control = control
    )

    if (!is.na(prev_exit) && length(ordered_tracts) > 1) {
      first <- ordered_tracts[1]
      last <- ordered_tracts[length(ordered_tracts)]
      d_first <- dist_val[paste(prev_exit, first, sep = "|")]
      d_last <- dist_val[paste(prev_exit, last, sep = "|")]
      if (!is.na(d_first) && !is.na(d_last) && d_last < d_first) {
        ordered_tracts <- rev(ordered_tracts)
      }
    }

    prev_exit <- ordered_tracts[length(ordered_tracts)]

    sequences[[j]] <- tibble::tibble(
      zone_id = zid,
      tract_id = ordered_tracts,
      visit_order = seq_along(ordered_tracts)
    )
  }

  dplyr::bind_rows(sequences)
}


#' Group Externally-Defined Zones Within a Distance Bound and Label Them
#'
#' Given zones that already exist (not built by [surveyzones_build_zones()]),
#' finds a center for each zone, groups zones within `D_max` of a group
#' center by reusing the same p-median machinery one level up, and assigns
#' each zone a stable `ggg.zzz` label from two independent partition-wide
#' TSPs. Also produces a tract-level visit sequence oriented across zone
#' boundaries.
#'
#' @param zone_assignments A data.frame with `tract_id`, `zone_id`, and
#'   optionally `partition_id`. Every tract in a given `zone_id` must
#'   share the same `partition_id`.
#' @param sparse_distances Tract-level sparse distance tibble (`origin_id`,
#'   `destination_id`, `distance`).
#' @param D_max Numeric scalar. Every zone in a group is within `D_max` of
#'   that group's center zone.
#' @param solver Character scalar. One of `"glpk"` (default), `"highs"`,
#'   `"cbc"`, `"symphony"`, `"spopt"` -- reused for the zone-level p-median.
#' @param max_time,rel_tol,verbose,use_cache Passed through to the reused
#'   solver, same meaning as in [surveyzones_build_zones()].
#' @param access_points An sf object with POINT geometries and a
#'   `tract_id` column. Used to fill missing pairwise distances via
#'   haversine. Required when `complete_distances = TRUE`.
#' @param speed_kmh Numeric scalar. Assumed travel speed for haversine
#'   fill-in (km/h). Default `0.1`.
#' @param complete_distances Logical scalar. When `TRUE` (default),
#'   missing pairwise distances are filled using haversine via
#'   [surveyzones_complete_distances()].
#' @param starting_points An optional sf POINT object with a
#'   `partition_id` column, one row per partition. When a partition has a
#'   matching row, that partition's `zzz` and `ggg` orders are anchored so
#'   the first zone and first group are the ones closest to it (see
#'   [.anchor_partition_order()]). Partitions absent from `starting_points`
#'   keep the unanchored TSP orientation. Requires `access_points`
#'   regardless of `complete_distances`.
#' @param method,control Seriation method/control passed to
#'   [surveyzones_sequence_tracts()] for tract-level ordering. Zone- and
#'   group-level TSPs always use `"TSP"`.
#'
#' @return A `surveyzones_zone_groups` object.
#'
#' @export
surveyzones_group_zones <- function(
  zone_assignments,
  sparse_distances,
  D_max,
  solver = "glpk",
  max_time = 300,
  rel_tol = 0.01,
  verbose = FALSE,
  use_cache = TRUE,
  access_points = NULL,
  speed_kmh = 0.1,
  complete_distances = TRUE,
  starting_points = NULL,
  method = "TSP",
  control = NULL
) {
  D_max <- as.numeric(D_max)
  validate_solver(solver)

  za <- .validate_zone_assignments(zone_assignments)

  if (!is.null(starting_points) && is.null(access_points)) {
    cli::cli_abort(c(
      "{.arg access_points} is required when {.arg starting_points} is supplied.",
      "i" = "Starting-point anchoring needs real coordinates to compute distances."
    ))
  }
  if (!is.null(starting_points) && !"partition_id" %in% names(starting_points)) {
    cli::cli_abort("{.arg starting_points} must contain a {.field partition_id} column.")
  }

  if (complete_distances) {
    if (is.null(access_points)) {
      cli::cli_abort(c(
        "{.arg access_points} is required when {.arg complete_distances} is {.val TRUE}.",
        "i" = "Either pass it explicitly or set {.code complete_distances = FALSE}."
      ))
    }
    sparse_distances <- surveyzones_complete_distances(
      sparse_distances, access_points, speed_kmh = speed_kmh
    )
  }

  centers <- .find_zone_centers(za, sparse_distances)

  partitions <- split(centers$zone_id, centers$partition_id)

  zones_out <- vector("list", length(partitions))
  groups_out <- vector("list", length(partitions))
  sequence_out <- vector("list", length(partitions))
  diagnostics_out <- vector("list", length(partitions))

  for (idx in seq_along(partitions)) {
    pid <- names(partitions)[idx]
    zone_ids <- partitions[[idx]]
    za_part <- za |> dplyr::filter(.data$partition_id == pid)
    centers_part <- centers |> dplyr::filter(.data$partition_id == pid)
    zone_centers_part <- stats::setNames(
      centers_part$center_tract_id, centers_part$zone_id
    )

    starting_point_part <- if (is.null(starting_points)) {
      NULL
    } else {
      sp <- starting_points[starting_points$partition_id == pid, ]
      if (nrow(sp) == 0) NULL else sp[1, ]
    }

    grouped <- .group_zones_pmedian(
      zone_ids = zone_ids, zone_assignments = za_part,
      sparse_distances = sparse_distances, D_max = D_max,
      solver = solver, max_time = max_time, rel_tol = rel_tol,
      verbose = verbose, use_cache = use_cache
    )

    labels <- .label_zone_groups(
      zone_ids = zone_ids, groups = grouped$groups,
      zone_assignments = za_part, sparse_distances = sparse_distances,
      zone_centers = zone_centers_part, starting_point = starting_point_part,
      access_points = access_points
    )

    zone_order <- labels$zone_labels |>
      dplyr::arrange(as.integer(.data$zzz)) |>
      dplyr::pull(.data$zone_id)

    sequence_out[[idx]] <- .sequence_group_tracts(
      zone_order, za_part, sparse_distances, method = method, control = control
    )

    part_zones <- centers |>
      dplyr::filter(.data$partition_id == pid) |>
      dplyr::left_join(grouped$assignments, by = "zone_id") |>
      dplyr::left_join(labels$zone_labels, by = "zone_id") |>
      dplyr::left_join(labels$group_labels, by = "group_id") |>
      dplyr::mutate(label = paste0(.data$ggg, ".", .data$zzz))

    zones_out[[idx]] <- part_zones |>
      dplyr::select(
        "zone_id", "partition_id", "center_tract_id", "n_tracts",
        "group_id", "ggg", "zzz", "label"
      )

    groups_out[[idx]] <- grouped$groups |>
      dplyr::left_join(labels$group_labels, by = "group_id") |>
      dplyr::mutate(partition_id = pid) |>
      dplyr::select(
        "group_id", "partition_id", "center_zone_id", "n_zones",
        "diameter", "ggg"
      )

    diagnostics_out[[idx]] <- grouped$diagnostics
  }

  new_surveyzones_zone_groups(
    zones = dplyr::bind_rows(zones_out),
    groups = dplyr::bind_rows(groups_out),
    sequence = dplyr::bind_rows(sequence_out),
    parameters = list(D_max = D_max, solver = solver, timestamp = Sys.time()),
    diagnostics = list(
      solver_status = purrr::map_chr(diagnostics_out, "solver_status"),
      objective_value = purrr::map_dbl(diagnostics_out, "objective_value"),
      n_variables = purrr::map_int(diagnostics_out, "n_variables"),
      solve_time = purrr::map_dbl(diagnostics_out, "solve_time"),
      optimality_gap = purrr::map_dbl(diagnostics_out, "optimality_gap")
    )
  )
}

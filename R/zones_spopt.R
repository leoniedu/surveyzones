# Internal spopt backend for surveyzones_build_zones
#
# Delegates zone-building to spopt::allocate_zones() and converts
# the result to the surveyzones assignments/zones format.


#' Solve a single partition using the spopt backend
#'
#' Dispatches to spopt::allocate_zones() using the uncap-then-split
#' or direct strategy, and returns the result in surveyzones format.
#'
#' @param full_sparse_distances Tibble of all distance pairs (unfiltered).
#' @param tracts Data frame with `tract_id`, `expected_service_time`.
#' @param D_max Numeric. Maximum distance.
#' @param max_workload_per_zone Numeric. Workload cap per zone.
#' @param strategy Character. Effective strategy: `"uncap_then_split"` or
#'   `"direct"`.
#' @param verbose Logical. Print progress?
#' @return A list with `assignments`, `zones`, and `diagnostics`.
#' @keywords internal
.solve_partition_spopt <- function(
  full_sparse_distances,
  tracts,
  D_max,
  max_workload_per_zone,
  strategy = "direct",
  verbose = FALSE
) {
  start_time <- proc.time()[["elapsed"]]
  capacitated <- is.finite(max_workload_per_zone)

  if (strategy == "uncap_then_split" && capacitated) {
    return(.spopt_uncap_then_split(
      full_sparse_distances = full_sparse_distances,
      tracts = tracts,
      D_max = D_max,
      max_workload_per_zone = max_workload_per_zone,
      verbose = verbose
    ))
  }

  # Direct solve: uncapacitated or capacitated
  tracts_sf <- .tracts_to_sf(tracts)

  method <- if (capacitated) "cflp" else "p_median"
  capacity <- if (capacitated) max_workload_per_zone else NULL

  cli::cli_alert_info(
    "spopt: method = {.val {method}}, D_max = {D_max}"
  )

  spopt_result <- tryCatch(
    spopt::allocate_zones(
      zones = tracts_sf,
      max_distance = D_max,
      method = method,
      weight_col = "expected_service_time",
      capacity = capacity,
      distances = full_sparse_distances,
      id_col = "tract_id",
      sequence = FALSE,
      verbose = verbose
    ),
    error = function(e) {
      cli::cli_alert_danger("spopt solve failed: {conditionMessage(e)}")
      NULL
    }
  )

  if (is.null(spopt_result)) {
    return(.empty_zone_solution(
      solver_status = "infeasible",
      n_variables = NA_integer_,
      solve_time = proc.time()[["elapsed"]] - start_time
    ))
  }

  elapsed <- proc.time()[["elapsed"]] - start_time
  result <- .spopt_result_to_zones(spopt_result, full_sparse_distances)

  cli::cli_alert_success(
    "spopt: {nrow(result$zones)} zones in {round(elapsed, 1)}s"
  )

  result$diagnostics <- list(
    solver_status = "optimal",
    objective_value = sum(result$assignments$distance_to_center, na.rm = TRUE),
    n_variables = NA_integer_,
    solve_time = elapsed,
    optimality_gap = NA_real_
  )

  result
}


#' Uncap-then-split strategy using spopt
#'
#' Phase 1: Solve uncapacitated p_median via spopt to get groups.
#' Phase 2: Split any oversized groups using spopt's cflp method.
#'
#' @inheritParams .solve_partition_spopt
#' @return A list with `assignments`, `zones`, and `diagnostics`.
#' @keywords internal
.spopt_uncap_then_split <- function(
  full_sparse_distances,
  tracts,
  D_max,
  max_workload_per_zone,
  verbose = FALSE
) {
  start_time <- proc.time()[["elapsed"]]
  tracts_sf <- .tracts_to_sf(tracts)

  n <- nrow(tracts)
  cli::cli_alert_info(
    "spopt uncap_then_split: Phase 1 \u2014 uncapacitated p_median ({n} tracts)"
  )

  # Phase 1: uncapacitated solve
  phase1_t0 <- proc.time()[["elapsed"]]
  uncap_result <- tryCatch(
    spopt::allocate_zones(
      zones = tracts_sf,
      max_distance = D_max,
      method = "p_median",
      weight_col = "expected_service_time",
      distances = full_sparse_distances,
      id_col = "tract_id",
      sequence = FALSE,
      verbose = verbose
    ),
    error = function(e) {
      cli::cli_alert_danger("spopt uncap solve failed: {conditionMessage(e)}")
      NULL
    }
  )

  if (is.null(uncap_result)) {
    return(.empty_zone_solution(
      solver_status = "infeasible",
      n_variables = NA_integer_,
      solve_time = proc.time()[["elapsed"]] - start_time
    ))
  }

  phase1_elapsed <- proc.time()[["elapsed"]] - phase1_t0
  n_groups <- nrow(uncap_result$zones |> dplyr::filter(.data$.is_center))
  cli::cli_alert_success(
    "Phase 1 done: {n_groups} groups in {round(phase1_elapsed, 1)}s"
  )

  uncap_zones <- .spopt_result_to_zones(uncap_result, full_sparse_distances)

  # Check workloads
  workload_by_tract <- stats::setNames(
    tracts$expected_service_time, as.character(tracts$tract_id)
  )

  zone_wl <- uncap_zones$assignments |>
    dplyr::mutate(.wl = workload_by_tract[.data$tract_id]) |>
    dplyr::summarise(total_wl = sum(.data$.wl), .by = "zone_id")

  oversized <- zone_wl |>
    dplyr::filter(.data$total_wl > max_workload_per_zone) |>
    dplyr::pull(.data$zone_id)

  if (length(oversized) == 0L) {
    cli::cli_alert_success("All zones satisfy the workload cap \u2014 no splitting needed")
    uncap_zones$assignments$group_id <- uncap_zones$assignments$zone_id
    uncap_zones$zones$group_id <- uncap_zones$zones$zone_id
    elapsed <- proc.time()[["elapsed"]] - start_time
    uncap_zones$diagnostics <- list(
      solver_status = "optimal",
      objective_value = sum(uncap_zones$assignments$distance_to_center, na.rm = TRUE),
      n_variables = NA_integer_,
      solve_time = elapsed,
      optimality_gap = NA_real_
    )
    return(uncap_zones)
  }

  n_good <- nrow(uncap_zones$zones) - length(oversized)
  cli::cli_alert_info(
    "Phase 2 \u2014 {n_good} zones OK, {length(oversized)} oversized zone{?s} to split"
  )

  # Partition into good and oversized
  good_assignments <- uncap_zones$assignments |>
    dplyr::filter(!.data$zone_id %in% oversized) |>
    dplyr::mutate(group_id = .data$zone_id)
  good_zones <- uncap_zones$zones |>
    dplyr::filter(!.data$zone_id %in% oversized) |>
    dplyr::mutate(group_id = .data$zone_id)

  split_assignments <- vector("list", length(oversized))
  split_zones <- vector("list", length(oversized))

  for (idx in seq_along(oversized)) {
    zid <- oversized[[idx]]
    zone_tract_ids <- uncap_zones$assignments$tract_id[
      uncap_zones$assignments$zone_id == zid
    ]
    zone_tracts <- tracts |> dplyr::filter(.data$tract_id %in% zone_tract_ids)
    zone_sf <- .tracts_to_sf(zone_tracts)

    zone_wl_total <- sum(workload_by_tract[zone_tract_ids])
    n_tracts_in_zone <- length(zone_tract_ids)

    cli::cli_alert_info(
      "  [{idx}/{length(oversized)}] Splitting zone {zid} ({n_tracts_in_zone} tracts, wl={round(zone_wl_total, 1)})"
    )

    sub_dist <- full_sparse_distances |>
      dplyr::filter(
        .data$origin_id %in% zone_tract_ids,
        .data$destination_id %in% zone_tract_ids
      )

    sub_result <- tryCatch(
      spopt::allocate_zones(
        zones = zone_sf,
        max_distance = D_max,
        method = "cflp",
        weight_col = "expected_service_time",
        capacity = max_workload_per_zone,
        distances = sub_dist,
        id_col = "tract_id",
        sequence = FALSE,
        verbose = verbose
      ),
      error = function(e) {
        cli::cli_alert_danger("  Split of zone {zid} failed: {conditionMessage(e)}")
        NULL
      }
    )

    if (is.null(sub_result)) {
      return(.empty_zone_solution(
        solver_status = "infeasible",
        n_variables = NA_integer_,
        solve_time = proc.time()[["elapsed"]] - start_time
      ))
    }

    sub_zones <- .spopt_result_to_zones(sub_result, full_sparse_distances)
    split_assignments[[idx]] <- sub_zones$assignments |>
      dplyr::mutate(group_id = zid)
    split_zones[[idx]] <- sub_zones$zones |>
      dplyr::mutate(group_id = zid)
  }

  final_assignments <- dplyr::bind_rows(
    good_assignments,
    dplyr::bind_rows(split_assignments)
  )
  final_zones <- dplyr::bind_rows(
    good_zones,
    dplyr::bind_rows(split_zones)
  )

  elapsed <- proc.time()[["elapsed"]] - start_time
  cli::cli_alert_success(
    "spopt uncap_then_split: {nrow(final_zones)} zones in {round(elapsed, 1)}s"
  )

  list(
    assignments = final_assignments,
    zones = final_zones,
    diagnostics = list(
      solver_status = "optimal",
      objective_value = sum(final_assignments$distance_to_center, na.rm = TRUE),
      n_variables = NA_integer_,
      solve_time = elapsed,
      optimality_gap = NA_real_
    )
  )
}


#' Convert spopt result to surveyzones assignments/zones format
#'
#' @param spopt_result The list returned by `spopt::allocate_zones()`.
#' @param full_sparse_distances Tibble of all distance pairs for diameter
#'   computation.
#' @return A list with `assignments` and `zones` tibbles.
#' @keywords internal
.spopt_result_to_zones <- function(spopt_result, full_sparse_distances) {
  zones_sf <- spopt_result$zones
  tract_ids <- as.character(zones_sf$tract_id)

  # Map center index to tract_id
  center_ids <- tract_ids[zones_sf$.center]

  assignments <- tibble::tibble(
    tract_id = tract_ids,
    zone_id = center_ids,
    partition_id = NA_character_,
    center_id = center_ids,
    distance_to_center = zones_sf$.distance
  )

  # Build zone summary
  workload <- stats::setNames(
    zones_sf$expected_service_time, tract_ids
  )

  zones <- assignments |>
    dplyr::summarise(
      center_tract_id = unique(.data$center_id),
      total_workload = sum(workload[.data$tract_id]),
      n_tracts = dplyr::n(),
      partition_id = NA_character_,
      .by = "zone_id"
    )

  zones$diameter <- purrr::map_dbl(
    zones$zone_id,
    \(zid) .compute_zone_diameter(
      assignments$tract_id[assignments$zone_id == zid],
      full_sparse_distances
    )
  )

  list(
    assignments = assignments |>
      dplyr::mutate(group_id = .data$zone_id) |>
      dplyr::select(
        "tract_id", "zone_id", "partition_id",
        "center_id", "distance_to_center", "group_id"
      ),
    zones = zones |>
      dplyr::mutate(group_id = .data$zone_id) |>
      dplyr::select(
        "zone_id", "partition_id", "center_tract_id",
        "total_workload", "diameter", "n_tracts", "group_id"
      )
  )
}


#' Create a minimal sf object from tracts
#'
#' spopt::allocate_zones() requires an sf object, but when pre-computed
#' distances are provided, geometry is never used. This helper creates
#' a dummy sf with POINT(0,0) geometries.
#'
#' @param tracts Data frame with `tract_id` and `expected_service_time`.
#' @return An sf object with point geometries.
#' @keywords internal
.tracts_to_sf <- function(tracts) {
  sf::st_as_sf(
    tracts,
    geometry = sf::st_sfc(
      rep(list(sf::st_point(c(0, 0))), nrow(tracts)),
      crs = 4326
    )
  )
}

#' Build Zones for Survey Tracts
#'
#' Top-level entry point that partitions tracts (optionally by
#' jurisdiction), solves a p-median problem for each partition,
#' and returns a combined `surveyzones_plan`.
#'
#' When `max_workload_per_zone` is finite, solves a **capacitated**
#' p-median (zones cannot exceed the workload cap).
#' When `max_workload_per_zone = Inf` (the default), solves an
#' **uncapacitated** p-median (no workload constraints, fewer MILP
#' constraints, tighter LP relaxation, faster solve times).
#'
#' @param sparse_distances A `data.table` of sparse distances as
#'   returned by [surveyzones_compute_sparse_distances()] or
#'   [surveyzones_precomputed_distances()].
#' @param tracts A data.frame with columns `tract_id` and
#'   `expected_service_time`.  Optionally `partition_id`.
#' @param D_max Numeric scalar.  Maximum distance.  Pairs in
#'   `sparse_distances` with `travel_time > D_max` are dropped
#'   before solving.
#' @param max_workload_per_zone Numeric scalar.  Upper bound on
#'   the sum of `expected_service_time` within any zone.  Defaults
#'   to `Inf` (uncapacitated).
#' @param target_zone_size Numeric scalar or `NULL`.  Desired number
#'   of tracts per zone, used to compute the initial number of zones
#'   K as `ceiling(n_tracts / target_zone_size)`.  When both
#'   `target_zone_size` and a finite `max_workload_per_zone` are
#'   given, K is the maximum of the two implied values.
#' @param K_max Integer.  Safety cap on the number of zones per
#'   partition.
#' @param enforce_partition Logical.  Whether to split by
#'   `partition_id` (default `TRUE`).
#' @param candidates Optional character vector of tract_ids eligible
#'   to be zone centers.  `NULL` (default) means all tracts are
#'   candidates.
#' @param max_variables Integer.  Safety limit on the number of
#'   x variables in the MILP.  Aborts if exceeded.
#' @param solver Character scalar.  Which MILP solver backend to use.
#'   One of `"glpk"` (default), `"highs"`, or `"cbc"`.  The
#'   corresponding ROI plugin package must be installed (e.g.,
#'   `ROI.plugin.highs`).
#'
#' @return A `surveyzones_plan` object.
#'
#' @export
surveyzones_build_zones <- function(
  sparse_distances,
  tracts,
  D_max,
  max_workload_per_zone = Inf,
  target_zone_size = NULL,
  K_max = NULL,
  enforce_partition = TRUE,
  candidates = NULL,
  max_variables = 500000L,
  solver = "glpk",
  max_time = 300,
  rel_tol = 0.01
) {
  validate_tracts(tracts)
  validate_solver(solver)

  if (is.null(target_zone_size) && is.infinite(max_workload_per_zone)) {
    cli::cli_abort(c(
      "Either {.arg target_zone_size} or a finite {.arg max_workload_per_zone} must be provided.",
      "i" = "Without at least one, the number of zones K cannot be determined."
    ))
  }

  # Keep unfiltered distances for true diameter computation
  full_sparse_distances <- sparse_distances

  # Filter distances to D_max
  n_before <- nrow(sparse_distances)
  sparse_distances <- sparse_distances[travel_time <= D_max]
  n_after <- nrow(sparse_distances)

  cli::cli_alert_info(
    "Input: {nrow(tracts)} tracts, {n_after} distance pairs (filtered from {n_before} by D_max = {D_max})"
  )
  workload_label <- if (is.infinite(max_workload_per_zone)) "Inf (uncapacitated)" else max_workload_per_zone
  size_label <- target_zone_size %||% "NULL"
  cli::cli_alert_info(
    "Parameters: D_max = {D_max}, max_workload = {workload_label}, target_zone_size = {size_label}, solver = {solver}"
  )

  parts <- surveyzones_partition(tracts, enforce_partition)
  cli::cli_alert_info(
    "{length(parts)} partition{?s}: {paste(names(parts), collapse = ', ')}"
  )

  all_assignments <- vector("list", length(parts))
  all_zones <- vector("list", length(parts))
  all_diag <- vector("list", length(parts))
  partition_names <- names(parts)

  build_t0 <- proc.time()[["elapsed"]]

  for (i in seq_along(parts)) {
    pid <- partition_names[i]
    part_tracts <- parts[[i]]

    cli::cli_h2(
      "Partition {i}/{length(parts)}: {pid} ({nrow(part_tracts)} tracts)"
    )

    # Filter distances to this partition's tracts
    part_ids <- as.character(part_tracts$tract_id)
    part_dist <- sparse_distances[
      origin_id %in% part_ids & destination_id %in% part_ids
    ]
    cli::cli_alert_info("Distance pairs for partition: {nrow(part_dist)}")

    part_candidates <- if (!is.null(candidates)) {
      intersect(candidates, part_ids)
    } else {
      part_ids
    }
    cli::cli_alert_info("Candidate centers: {length(part_candidates)}")

    k_max_part <- K_max %||% nrow(part_tracts)

    # Full (unfiltered) distances for this partition, used for diameter
    part_full_dist <- full_sparse_distances[
      origin_id %in% part_ids & destination_id %in% part_ids
    ]

    result <- surveyzones_build_zones_single(
      sparse_distances = part_dist,
      full_sparse_distances = part_full_dist,
      tracts = part_tracts,
      D_max = D_max,
      max_workload_per_zone = max_workload_per_zone,
      target_zone_size = target_zone_size,
      K_max = k_max_part,
      candidates = part_candidates,
      max_variables = max_variables,
      solver = solver,
      max_time = max_time,
      rel_tol = rel_tol
    )

    # Tag with partition_id
    result$assignments$partition_id <- pid
    result$zones$partition_id <- pid

    all_assignments[[i]] <- result$assignments
    all_zones[[i]] <- result$zones
    all_diag[[i]] <- result$diagnostics

    cli::cli_alert_info(
      "Partition {pid}: status={result$diagnostics$solver_status}, zones={nrow(result$zones)}, time={round(result$diagnostics$solve_time, 2)}s"
    )
  }

  build_elapsed <- proc.time()[["elapsed"]] - build_t0

  assignments <- do.call(rbind, all_assignments)
  zones <- do.call(rbind, all_zones)

  cli::cli_rule()
  cli::cli_alert_success(
    "All partitions done: {nrow(zones)} zones for {nrow(assignments)} tracts in {round(build_elapsed, 2)}s"
  )
  diagnostics <- list(
    solver_status = vapply(all_diag, \(d) d$solver_status, character(1)),
    objective_value = vapply(all_diag, \(d) d$objective_value, numeric(1)),
    n_variables = vapply(all_diag, \(d) d$n_variables, integer(1)),
    solve_time = vapply(all_diag, \(d) d$solve_time, numeric(1))
  )

  parameters <- list(
    D_max = D_max,
    max_workload_per_zone = max_workload_per_zone,
    target_zone_size = target_zone_size,
    K_max = K_max,
    enforce_partition = enforce_partition,
    solver = solver,
    timestamp = Sys.time()
  )

  new_surveyzones_plan(assignments, zones, parameters, diagnostics)
}


#' Solve Zones for a Single Partition
#'
#' Tries increasing values of K (number of zones) until a feasible
#' solution is found.
#'
#' @inheritParams surveyzones_build_zones
#' @param candidates Character vector of candidate center tract_ids.
#'
#' @return A list with `assignments`, `zones`, and `diagnostics`.
#' @keywords internal
surveyzones_build_zones_single <- function(
  sparse_distances,
  full_sparse_distances,
  tracts,
  D_max,
  max_workload_per_zone = Inf,
  target_zone_size = NULL,
  K_max,
  candidates,
  max_variables = 500000L,
  solver = "glpk",
  max_time = 300,
  rel_tol = 0.01
) {
  n <- nrow(tracts)
  total_workload <- sum(tracts$expected_service_time)

  # Compute initial K from target_zone_size and/or workload cap
  K_from_size <- if (!is.null(target_zone_size)) {
    ceiling(n / target_zone_size)
  } else {
    1L
  }

  K_from_workload <- if (is.finite(max_workload_per_zone)) {
    ceiling(total_workload / max_workload_per_zone)
  } else {
    1L
  }

  K0 <- max(1L, K_from_size, K_from_workload)

  capacitated <- is.finite(max_workload_per_zone)
  model_type <- if (capacitated) "capacitated" else "uncapacitated"

  cli::cli_alert_info(
    "n = {n}, total workload = {round(total_workload, 1)}, K0 = {K0}, K_max = {K_max}, model = {model_type}"
  )

  for (K in seq(K0, K_max)) {
    cli::cli_alert("Trying K = {K}...")

    result <- surveyzones_solve_fixed_K(
      sparse_distances = sparse_distances,
      full_sparse_distances = full_sparse_distances,
      tracts = tracts,
      K = K,
      max_workload_per_zone = max_workload_per_zone,
      candidates = candidates,
      max_variables = max_variables,
      solver = solver,
      max_time = max_time,
      rel_tol = rel_tol
    )

    if (result$diagnostics$solver_status == "optimal") {
      cli::cli_alert_success(
        "Feasible at K = {K} (objective = {round(result$diagnostics$objective_value, 2)})"
      )
      return(result)
    }

    cli::cli_alert_warning("Infeasible at K = {K}, incrementing...")
  }

  cli::cli_alert_danger("No feasible solution found up to K_max = {K_max}.")
  # Return an infeasible result
  list(
    assignments = tibble::tibble(
      tract_id = character(0),
      zone_id = character(0),
      partition_id = character(0),
      center_id = character(0),
      travel_time_to_center = numeric(0)
    ),
    zones = tibble::tibble(
      zone_id = character(0),
      partition_id = character(0),
      center_tract_id = character(0),
      total_workload = numeric(0),
      diameter = numeric(0),
      n_tracts = integer(0)
    ),
    diagnostics = list(
      solver_status = "infeasible",
      objective_value = NA_real_,
      n_variables = NA_integer_,
      solve_time = NA_real_
    )
  )
}


#' Solve the Capacitated p-Median Problem for Fixed K
#'
#' Formulates and solves a mixed-integer linear program:
#' assign each tract to exactly one center, open exactly K centers,
#' respect workload capacity, and minimise total weighted distance.
#'
#' @inheritParams surveyzones_build_zones
#' @param K Integer.  Number of zones (centers) to open.
#' @param candidates Character vector of candidate center tract_ids.
#'
#' @return A list with `assignments`, `zones`, and `diagnostics`.
#' @keywords internal
surveyzones_solve_fixed_K <- function(
  sparse_distances,
  full_sparse_distances,
  tracts,
  K,
  max_workload_per_zone,
  candidates,
  max_variables = 500000L,
  solver = "glpk",
  max_time = 300,
  rel_tol = 0.01
) {
  start_time <- proc.time()["elapsed"]

  tract_ids <- as.character(tracts$tract_id)
  n <- length(tract_ids)

  cli::cli_alert_info("  K={K}: {n} tracts, building MILP...")

  # Build workload lookup
  workload <- stats::setNames(tracts$expected_service_time, tract_ids)

  # Candidate centers: indices into tract_ids
  cand_idx <- which(tract_ids %in% candidates)
  m <- length(cand_idx)

  if (m == 0L) {
    cli::cli_abort("No candidate centers available for this partition.")
  }

  # Sparse pairs: (tract index i, candidate index j) where distance exists
  dt <- data.table::copy(sparse_distances)
  # Map IDs to indices
  tract_to_idx <- stats::setNames(seq_len(n), tract_ids)
  cand_to_jdx <- stats::setNames(seq_along(cand_idx), tract_ids[cand_idx])

  # Keep only pairs where destination is a candidate
  dt <- dt[destination_id %in% tract_ids[cand_idx]]
  # Also allow self-assignment (center assigned to itself with distance 0)
  self_pairs <- data.table::data.table(
    origin_id = tract_ids[cand_idx],
    destination_id = tract_ids[cand_idx],
    travel_time = 0
  )
  dt <- data.table::rbindlist(list(dt, self_pairs))
  dt <- unique(dt, by = c("origin_id", "destination_id"))

  # Map to numeric indices
  dt[, i := tract_to_idx[origin_id]]
  dt[, j := cand_to_jdx[destination_id]]
  dt <- dt[!is.na(i) & !is.na(j)]

  n_x <- nrow(dt)

  if (n_x > max_variables) {
    cli::cli_abort(c(
      "Number of x-variables ({n_x}) exceeds {.arg max_variables} ({max_variables}).",
      "i" = "Reduce D_max or partition the problem further."
    ))
  }

  capacitated <- is.finite(max_workload_per_zone)
  n_cap <- if (capacitated) m else 0L

  cli::cli_alert_info(
    "  K={K}: {n_x} x-variables, {m} candidates, {n} assignment constraints, {n_cap} capacity constraints"
  )

  # Check that every tract can reach at least one candidate
  reachable <- unique(dt$i)
  unreachable <- setdiff(seq_len(n), reachable)
  if (length(unreachable) > 0) {
    unreach_ids <- tract_ids[unreachable]
    cli::cli_alert_warning(
      "{length(unreachable)} tract{?s} unreachable: {.val {head(unreach_ids, 5)}}{if (length(unreach_ids) > 5) '...'}"
    )
    return(list(
      assignments = tibble::tibble(
        tract_id = character(0),
        zone_id = character(0),
        partition_id = character(0),
        center_id = character(0),
        travel_time_to_center = numeric(0)
      ),
      zones = tibble::tibble(
        zone_id = character(0),
        partition_id = character(0),
        center_tract_id = character(0),
        total_workload = numeric(0),
        diameter = numeric(0),
        n_tracts = integer(0)
      ),
      diagnostics = list(
        solver_status = "infeasible",
        objective_value = NA_real_,
        n_variables = as.integer(n_x),
        solve_time = as.numeric(proc.time()["elapsed"] - start_time)
      )
    ))
  }

  # Build the MILP with ompr
  # Sparse x: only for pairs that exist in dt
  # Index the sparse pairs as 1..n_x
  pair_i <- dt$i
  pair_j <- dt$j
  pair_d <- dt$travel_time

  # Workload for each tract
  w <- workload[tract_ids]

  model <- ompr::MILPModel() |>
    ompr::add_variable(x[p], p = 1:n_x, type = "binary") |>
    ompr::add_variable(y[j], j = 1:m, type = "binary") |>
    ompr::set_objective(
      ompr::sum_over(pair_d[p] * x[p], p = 1:n_x),
      sense = "min"
    ) |>
    # x[p] <= y[j] — can only assign to an open center
    ompr::add_constraint(
      x[p] <= y[pair_j[p]],
      p = 1:n_x
    ) |>
    # Exactly K centers
    ompr::add_constraint(
      ompr::sum_over(y[j], j = 1:m) == K
    )

  # Each tract assigned exactly once (loop: ompr sum_over doesn't

  # support external vector indexing in quantified constraints)
  for (i_val in seq_len(n)) {
    ps <- which(pair_i == i_val)
    model <- model |>
      ompr::add_constraint(ompr::sum_over(x[p], p = ps) == 1)
  }

  # Workload capacity per center (skip for uncapacitated model)
  if (capacitated) {
    for (j_val in seq_len(m)) {
      ps <- which(pair_j == j_val)
      if (length(ps) > 0) {
        model <- model |>
          ompr::add_constraint(
            ompr::sum_over(w[pair_i[p]] * x[p], p = ps) <= max_workload_per_zone
          )
      }
    }
  }

  # Solve — solver-specific parameter dispatch (cf. orce package)
  solver_label <- toupper(solver)
  cli::cli_alert_info(
    "  K={K}: Solving MILP with {solver_label} (max_time={max_time}s, rel_tol={rel_tol})..."
  )
  solve_t0 <- proc.time()[["elapsed"]]
  if (solver == "symphony") {
    result <- ompr::solve_model(
      model,
      ompr.roi::with_ROI(
        solver = solver,
        max_time = as.numeric(max_time),
        gap_limit = rel_tol * 100
      )
    )
  } else {
    result <- ompr::solve_model(
      model,
      ompr.roi::with_ROI(
        solver = solver,
        max_time = as.numeric(max_time),
        rel_tol = rel_tol
      )
    )
  }
  solve_elapsed <- proc.time()[["elapsed"]] - solve_t0

  solve_time <- as.numeric(proc.time()["elapsed"] - start_time)
  status <- ompr::solver_status(result)
  # ROI returns "success" while other solvers return "optimal"
  if (status == "success") {
    status <- "optimal"
  }
  # Accept feasible solutions found within time limit
  if (status == "error") {
    has_solution <- tryCatch(
      { ompr::objective_value(result); TRUE },
      error = function(e) FALSE
    )
    if (has_solution) {
      cli::cli_alert_warning(
        "  K={K}: solver hit time/gap limit but found a feasible solution"
      )
      status <- "optimal"
    }
  }

  cli::cli_alert_info(
    "  K={K}: solver returned {.val {status}} in {round(solve_elapsed, 2)}s (total prep+solve: {round(solve_time, 2)}s)"
  )

  if (status != "optimal") {
    return(list(
      assignments = tibble::tibble(
        tract_id = character(0),
        zone_id = character(0),
        partition_id = character(0),
        center_id = character(0),
        travel_time_to_center = numeric(0)
      ),
      zones = tibble::tibble(
        zone_id = character(0),
        partition_id = character(0),
        center_tract_id = character(0),
        total_workload = numeric(0),
        diameter = numeric(0),
        n_tracts = integer(0)
      ),
      diagnostics = list(
        solver_status = status,
        objective_value = NA_real_,
        n_variables = as.integer(n_x),
        solve_time = solve_time
      )
    ))
  }

  # Extract solution
  x_sol <- ompr::get_solution(result, x[p])
  active <- x_sol[x_sol$value > 0.5, ]

  cli::cli_alert_info(
    "  K={K}: {nrow(active)} active assignments, objective = {round(ompr::objective_value(result), 2)}"
  )

  assignments <- tibble::tibble(
    tract_id = tract_ids[pair_i[active$p]],
    center_id = tract_ids[cand_idx[pair_j[active$p]]],
    travel_time_to_center = pair_d[active$p]
  )
  assignments$zone_id <- assignments$center_id
  assignments$partition_id <- NA_character_

  # Zone summary
  zones <- assignments |>
    dplyr::summarise(
      center_tract_id = unique(center_id),
      total_workload = sum(workload[tract_id]),
      n_tracts = dplyr::n(),
      .by = zone_id
    )
  zones$partition_id <- NA_character_

  # Diameter: max pairwise distance among tracts in each zone.
  # Uses full (unfiltered) sparse distances so peripheral pairs beyond
  # D_max are included.  Returns NA when any member pair is missing
  # from the distance table (true diameter unknown).
  zone_members <- split(assignments$tract_id, assignments$zone_id)
  zones$diameter <- vapply(zones$zone_id, function(zid) {
    members <- zone_members[[zid]]
    n_members <- length(members)
    if (n_members <= 1L) return(0)
    intra <- full_sparse_distances[
      origin_id %in% members & destination_id %in% members &
        origin_id != destination_id
    ]
    if (nrow(intra) == 0L) return(NA_real_)
    # Check completeness: need all n*(n-1)/2 unique unordered pairs
    pairs <- unique(paste0(
      pmin(intra$origin_id, intra$destination_id), "|",
      pmax(intra$origin_id, intra$destination_id)
    ))
    expected <- n_members * (n_members - 1L) %/% 2L
    if (length(pairs) < expected) return(NA_real_)
    max(intra$travel_time)
  }, numeric(1))

  cli::cli_alert_info(
    "  K={K}: {nrow(zones)} zones, workload range [{round(min(zones$total_workload), 1)}-{round(max(zones$total_workload), 1)}], max diameter {round(max(zones$diameter, na.rm = TRUE), 2)} km"
  )

  list(
    assignments = assignments[, c(
      "tract_id",
      "zone_id",
      "partition_id",
      "center_id",
      "travel_time_to_center"
    )],
    zones = zones[, c(
      "zone_id",
      "partition_id",
      "center_tract_id",
      "total_workload",
      "diameter",
      "n_tracts"
    )],
    diagnostics = list(
      solver_status = status,
      objective_value = ompr::objective_value(result),
      n_variables = as.integer(n_x),
      solve_time = solve_time
    )
  )
}


#' Validate and Load a Solver Plugin
#'
#' Checks that the requested solver is supported and that the
#' corresponding ROI plugin package is available.
#'
#' @param solver Character scalar: `"glpk"`, `"highs"`, or `"cbc"`.
#' @param call Caller environment for error reporting.
#' @return `solver`, invisibly.
#' @keywords internal
validate_solver <- function(solver, call = rlang::caller_env()) {
  supported <- c("glpk", "highs", "cbc", "symphony")
  if (!is.character(solver) || length(solver) != 1L) {
    cli::cli_abort(
      "{.arg solver} must be a single character string.",
      call = call
    )
  }
  if (!solver %in% supported) {
    cli::cli_abort(
      c(
        "{.arg solver} must be one of {.val {supported}}, not {.val {solver}}.",
        "i" = "Install the matching ROI plugin: {.pkg ROI.plugin.{solver}}"
      ),
      call = call
    )
  }
  pkg <- paste0("ROI.plugin.", solver)
  if (!requireNamespace(pkg, quietly = TRUE)) {
    cli::cli_abort(
      c(
        "Solver {.val {solver}} requires package {.pkg {pkg}}.",
        "i" = "Install it with: {.code install.packages(\"{pkg}\")}"
      ),
      call = call
    )
  }
  invisible(solver)
}


#' Validate Tract Table
#'
#' @param tracts Object to validate.
#' @param call Caller environment for error reporting.
#' @return `tracts`, invisibly.
#' @keywords internal
validate_tracts <- function(tracts, call = rlang::caller_env()) {
  if (!is.data.frame(tracts)) {
    cli::cli_abort("{.arg tracts} must be a data.frame.", call = call)
  }

  required <- c("tract_id", "expected_service_time")
  missing_cols <- setdiff(required, names(tracts))
  if (length(missing_cols) > 0) {
    cli::cli_abort(
      "{.arg tracts} is missing column{?s}: {.val {missing_cols}}.",
      call = call
    )
  }

  if (anyDuplicated(tracts$tract_id)) {
    cli::cli_abort(
      "{.field tract_id} values in {.arg tracts} must be unique.",
      call = call
    )
  }

  if (!is.numeric(tracts$expected_service_time)) {
    cli::cli_abort(
      "{.field expected_service_time} must be numeric.",
      call = call
    )
  }

  if (any(tracts$expected_service_time <= 0, na.rm = TRUE)) {
    cli::cli_abort(
      "{.field expected_service_time} must be > 0 for all tracts.",
      call = call
    )
  }

  invisible(tracts)
}

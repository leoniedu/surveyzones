utils::globalVariables(c("center_tract_id", "total_workload", "diameter", "n_tracts"))

#' @keywords internal
.is_solution_status <- function(status) {
  status %in% c("optimal", "feasible")
}

#' @keywords internal
.normalize_solver_status <- function(result) {
  raw_status <- ompr::solver_status(result)
  status <- if (raw_status == "success") "optimal" else raw_status
  obj <- tryCatch(ompr::objective_value(result), error = function(e) NA_real_)
  if (!.is_solution_status(status) && is.finite(obj)) {
    status <- "feasible"
  }
  optimality_gap <- tryCatch(
    result$additional_solver_output$ROI$message$info$mip_gap,
    error = function(e) NA_real_
  )
  if (is.null(optimality_gap)) optimality_gap <- NA_real_
  list(
    status = status,
    raw_status = raw_status,
    objective_value = obj,
    optimality_gap = as.numeric(optimality_gap)
  )
}

#' @keywords internal
.solve_model_with_status <- function(model, solver, max_time, rel_tol, verbose) {
  params <- .solver_control_params(solver, max_time, rel_tol, verbose)
  result <- ompr::solve_model(
    model,
    do.call(ompr.roi::with_ROI, params)
  )
  info <- .normalize_solver_status(result)
  c(list(result = result), info)
}

#' @keywords internal
.empty_zone_solution <- function(
  solver_status = "infeasible",
  objective_value = NA_real_,
  n_variables = NA_integer_,
  solve_time = NA_real_
) {
  list(
    assignments = tibble::tibble(
      tract_id = character(0),
      zone_id = character(0),
      partition_id = character(0),
      center_id = character(0),
      distance_to_center = numeric(0),
      group_id = character(0)
    ),
    zones = tibble::tibble(
      zone_id = character(0),
      partition_id = character(0),
      center_tract_id = character(0),
      total_workload = numeric(0),
      diameter = numeric(0),
      n_tracts = integer(0),
      group_id = character(0)
    ),
    diagnostics = list(
      solver_status = solver_status,
      objective_value = objective_value,
      n_variables = n_variables,
      solve_time = solve_time,
      optimality_gap = NA_real_
    )
  )
}

#' @keywords internal
.surveyzones_build_zones_impl <- function(
  sparse_distances,
  tracts,
  D_max,
  max_workload_per_zone = Inf,
  enforce_partition = TRUE,
  solver = "glpk",
  max_time = 300,
  rel_tol = 0.01,
  verbose = FALSE,
  strategy = "auto",
  use_cache = TRUE
) {
  # Coerce to numeric — integer D_max can cause solver infeasibility
  D_max <- as.numeric(D_max)
  max_workload_per_zone <- as.numeric(max_workload_per_zone)

  validate_tracts(tracts)
  validate_solver(solver)

  solve_single <- if (use_cache) {
    surveyzones_build_zones_single_mem
  } else {
    surveyzones_build_zones_single
  }

  # Count distances before and after D_max filtering
  n_before <- nrow(sparse_distances)
  n_after <- nrow(sparse_distances |> dplyr::filter(distance <= D_max))

  cli::cli_alert_info(
    "Input: {nrow(tracts)} tracts, {n_after} distance pairs (filtered from {n_before} by D_max = {D_max})"
  )
  workload_label <- if (is.infinite(max_workload_per_zone)) "Inf (uncapacitated)" else max_workload_per_zone
  cli::cli_alert_info(
    "Parameters: D_max = {D_max}, max_workload = {workload_label}, solver = {solver}"
  )

  parts <- surveyzones_partition(tracts, enforce_partition)
  cli::cli_alert_info(
    "{length(parts)} partition{?s}: {paste(names(parts), collapse = ', ')}"
  )

  build_t0 <- proc.time()[["elapsed"]]

  # Process each partition using purrr::imap
  results <- purrr::imap(parts, \(part_tracts, pid) {
    part_idx <- match(pid, names(parts))
    cli::cli_h2(
      "Partition {part_idx}/{length(parts)}: {pid} ({nrow(part_tracts)} tracts)"
    )

    # Unfiltered distances for this partition
    part_ids <- as.character(part_tracts$tract_id)
    part_full_dist <- sparse_distances |>
      dplyr::filter(origin_id %in% part_ids, destination_id %in% part_ids)
    cli::cli_alert_info("Distance pairs for partition (unfiltered): {nrow(part_full_dist)}")

    cli::cli_alert_info("Candidate centers: {length(part_ids)} (all tracts)")

    # Report cache status for this partition
    if (use_cache) {
      single_args <- list(
        full_sparse_distances = part_full_dist,
        tracts = part_tracts,
        D_max = D_max,
        max_workload_per_zone = max_workload_per_zone,
        solver = solver,
        max_time = max_time,
        rel_tol = rel_tol,
        verbose = verbose,
        strategy = strategy
      )
      is_cached <- do.call(
        memoise::has_cache(surveyzones_build_zones_single_mem),
        single_args
      )
      if (is_cached) {
        cli::cli_alert_success("Partition {pid}: using cached result.")
      } else {
        cli::cli_alert_info("Partition {pid}: computing and caching result.")
      }
    }

    result <- solve_single(
      full_sparse_distances = part_full_dist,
      tracts = part_tracts,
      D_max = D_max,
      max_workload_per_zone = max_workload_per_zone,
      solver = solver,
      max_time = max_time,
      rel_tol = rel_tol,
      verbose = verbose,
      strategy = strategy
    )

    # Tag with partition_id
    result$assignments$partition_id <- pid
    result$zones$partition_id <- pid

    cli::cli_alert_info(
      "Partition {pid}: status={result$diagnostics$solver_status}, zones={nrow(result$zones)}, time={round(result$diagnostics$solve_time, 2)}s"
    )

    result
  })

  build_elapsed <- proc.time()[["elapsed"]] - build_t0

  # Extract and combine results using dplyr::bind_rows
  assignments <- purrr::map(results, \(r) r$assignments) |> purrr::list_rbind()
  zones <- purrr::map(results, \(r) r$zones) |> purrr::list_rbind()

  cli::cli_rule()
  cli::cli_alert_success(
    "All partitions done: {nrow(zones)} zones for {nrow(assignments)} tracts in {round(build_elapsed, 2)}s"
  )

  # Extract diagnostics from each partition's results
  diagnostics <- list(
    solver_status = purrr::map_chr(results, \(r) r$diagnostics$solver_status),
    objective_value = purrr::map_dbl(results, \(r) r$diagnostics$objective_value),
    n_variables = purrr::map_int(results, \(r) r$diagnostics$n_variables),
    solve_time = purrr::map_dbl(results, \(r) r$diagnostics$solve_time),
    optimality_gap = purrr::map_dbl(results, \(r) r$diagnostics$optimality_gap)
  )

  parameters <- list(
    D_max = D_max,
    max_workload_per_zone = max_workload_per_zone,
    enforce_partition = enforce_partition,
    solver = solver,
    timestamp = Sys.time()
  )

  new_surveyzones_plan(assignments, zones, parameters, diagnostics)
}


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
#' Results are cached to disk by default (`use_cache = TRUE`): repeated calls
#' with identical arguments return instantly without re-running the solver.
#' Use `surveyzones_clear_cache()` to wipe the cache.
#'
#' @param sparse_distances A tibble of sparse distances as
#'   returned by [surveyzones_compute_sparse_distances()] or
#'   [surveyzones_precomputed_distances()].
#' @param tracts A data.frame with columns `tract_id` and
#'   `expected_service_time`.  Optionally `partition_id`.
#' @param D_max Numeric scalar.  Maximum distance.  Pairs in
#'   `sparse_distances` with `distance > D_max` are dropped
#'   before solving.
#' @param max_workload_per_zone Numeric scalar.  Upper bound on
#'   the sum of `expected_service_time` within any zone.  Defaults
#'   to `Inf` (uncapacitated).
#' @param enforce_partition Logical.  Whether to split by
#'   `partition_id` (default `TRUE`).
#' @param solver Character scalar.  Which solver backend to use.
#'   MILP solvers: `"glpk"` (default), `"highs"`, `"cbc"`, `"symphony"`
#'   (requires the corresponding ROI plugin, e.g. `ROI.plugin.highs`).
#'   Use `"spopt"` to delegate to [spopt::allocate_zones()] which uses
#'   a Rust/HiGHS backend (much faster for large problems; requires the
#'   `spopt` package).
#' @param max_time Numeric scalar.  Maximum seconds to spend solving each
#'   partition. Default `300`.
#' @param rel_tol Numeric scalar.  Relative tolerance for solver convergence.
#'   Default `0.01`.
#' @param verbose Logical.  Print solver progress messages? Default `FALSE`.
#' @param strategy Character scalar.  Solving strategy for capacitated problems
#'   (when `max_workload_per_zone` is finite).  `"auto"` (default) uses
#'   `"uncap_then_split"` when `max_workload_per_zone` is finite and `"direct"`
#'   otherwise.  `"uncap_then_split"` solves the uncapacitated p-median first
#'   (fast), then recursively splits any zone that exceeds the workload cap into
#'   smaller zones -- each split sub-problem is tiny and solved instantly.
#'   `"direct"` uses the capacitated MILP directly (may be slow for large
#'   instances).
#' @param use_cache Logical.  Whether to use the disk cache (default `TRUE`).
#'   Set to `FALSE` to force re-computation.
#' @param access_points Optional sf object with POINT geometries and a
#'   `tract_id` column.  Stored in the returned plan so that downstream
#'   functions like [surveyzones_sequence()] can use it without
#'   the caller having to pass it again.
#'
#' @return A `surveyzones_plan` object.
#'
#' @export
surveyzones_build_zones <- function(
  sparse_distances,
  tracts,
  D_max,
  max_workload_per_zone = Inf,
  enforce_partition = TRUE,
  solver = "glpk",
  max_time = 300,
  rel_tol = 0.01,
  verbose = FALSE,
  strategy = "auto",
  use_cache = TRUE,
  access_points = NULL
) {
  plan <- .surveyzones_build_zones_impl(
    sparse_distances      = sparse_distances,
    tracts                = tracts,
    D_max                 = D_max,
    max_workload_per_zone = max_workload_per_zone,
    enforce_partition     = enforce_partition,
    solver                = solver,
    max_time              = max_time,
    rel_tol               = rel_tol,
    verbose               = verbose,
    strategy              = strategy,
    use_cache             = use_cache
  )

  # Attach access points for downstream use (not part of cache key)
  plan$access_points <- access_points
  plan
}


#' Solve Zones for a Single Partition
#'
#' Tries increasing values of K (number of zones) until a feasible
#' solution is found.
#'
#' @inheritParams surveyzones_build_zones
#' @param full_sparse_distances Tibble of all distance pairs (unfiltered).
#'
#' @return A list with `assignments`, `zones`, and `diagnostics`.
#' @keywords internal
surveyzones_build_zones_single <- function(
  full_sparse_distances,
  tracts,
  D_max,
  max_workload_per_zone = Inf,
  solver = "glpk",
  max_time = 300,
  rel_tol = 0.01,
  verbose = FALSE,
  strategy = "auto"
) {
  n <- nrow(tracts)
  K_max <- n
  total_workload <- sum(tracts$expected_service_time)

  capacitated <- is.finite(max_workload_per_zone)

  valid_strategies <- c("auto", "uncap_then_split", "direct")
  if (!strategy %in% valid_strategies) {
    cli::cli_abort(
      "{.arg strategy} must be one of {.val {valid_strategies}}, not {.val {strategy}}."
    )
  }

  # Resolve effective strategy early so K0 can depend on it
  effective_strategy <- if (strategy == "auto") {
    if (capacitated) "uncap_then_split" else "direct"
  } else {
    strategy
  }

  # For uncap_then_split, start K0 = 1: the split phase handles capacity.
  # K_from_workload assumes even distribution and overestimates K when
  # population is clustered (e.g. bimodal).
  # For "direct", use K_from_workload to skip obviously infeasible K values.
  K0 <- if (effective_strategy == "direct" && capacitated) {
    max(1L, ceiling(total_workload / max_workload_per_zone))
  } else {
    1L
  }
  model_type <- if (capacitated) "capacitated" else "uncapacitated"

  cli::cli_alert_info(
    "n = {n}, total workload = {round(total_workload, 1)}, K0 = {K0}, model = {model_type}"
  )

  # -- Connected component decomposition --------------------------------------
  # Filter distances to this partition and D_max for graph connectivity
  tract_ids_chr <- as.character(tracts$tract_id)
  sparse_distances <- full_sparse_distances |>
    dplyr::filter(
      distance <= D_max,
      origin_id %in% tract_ids_chr,
      destination_id %in% tract_ids_chr
    )

  # If the sparse distance graph has multiple disconnected components,
  # solve each independently (much smaller MILPs).
  g <- igraph::graph_from_data_frame(
    data.frame(
      from = sparse_distances$origin_id,
      to = sparse_distances$destination_id,
      stringsAsFactors = FALSE
    ),
    directed = FALSE,
    vertices = data.frame(name = tract_ids_chr, stringsAsFactors = FALSE)
  )
  comp <- igraph::components(g)

  if (comp$no > 1) {
    comp_sizes <- sort(comp$csize, decreasing = TRUE)
    cli::cli_alert_info(
      "Disconnected graph: {comp$no} components (sizes: {paste(comp_sizes, collapse = ', ')}) \u2014 solving each independently"
    )

    # Map each tract to its component
    comp_membership <- stats::setNames(comp$membership, igraph::V(g)$name)

    all_comp_assignments <- vector("list", comp$no)
    all_comp_zones <- vector("list", comp$no)
    total_obj <- 0
    overall_comp_status <- "optimal"
    start_time <- proc.time()[["elapsed"]]

    for (ci in seq_len(comp$no)) {
      comp_ids <- names(comp_membership[comp_membership == ci])
      comp_n <- length(comp_ids)

      cli::cli_alert_info("Component {ci}/{comp$no}: {comp_n} tracts")

      comp_tracts <- tracts |> dplyr::filter(tract_id %in% comp_ids)
      comp_full_dist <- full_sparse_distances |>
        dplyr::filter(origin_id %in% comp_ids, destination_id %in% comp_ids)
      # Solve this component (will find 1 component and use normal K loop)
      comp_result <- surveyzones_build_zones_single(
        full_sparse_distances = comp_full_dist,
        tracts = comp_tracts,
        D_max = D_max,
        max_workload_per_zone = max_workload_per_zone,
        solver = solver,
        max_time = max_time,
        rel_tol = rel_tol,
        verbose = verbose,
        strategy = strategy
      )

      all_comp_assignments[[ci]] <- comp_result$assignments
      all_comp_zones[[ci]] <- comp_result$zones
      if (!.is_solution_status(comp_result$diagnostics$solver_status)) {
        bad_status <- comp_result$diagnostics$solver_status
        if (is.null(bad_status) || is.na(bad_status) || !nzchar(bad_status)) {
          bad_status <- "infeasible"
        }
        cli::cli_alert_danger(
          "Component {ci}/{comp$no} failed (status: {comp_result$diagnostics$solver_status}); partition infeasible"
        )
        return(.empty_zone_solution(
          solver_status = bad_status,
          n_variables = NA_integer_,
          solve_time = proc.time()[["elapsed"]] - start_time
        ))
      }
      if (comp_result$diagnostics$solver_status == "feasible") {
        overall_comp_status <- "feasible"
      }
      if (!is.na(comp_result$diagnostics$objective_value)) {
        total_obj <- total_obj + comp_result$diagnostics$objective_value
      }
    }

    combined_assignments <- dplyr::bind_rows(all_comp_assignments)
    combined_zones <- dplyr::bind_rows(all_comp_zones)
    elapsed <- proc.time()[["elapsed"]] - start_time

    cli::cli_alert_success(
      "All {comp$no} subgraphs solved: {nrow(combined_zones)} zones in {round(elapsed, 1)}s"
    )

    return(list(
      assignments = combined_assignments,
      zones = combined_zones,
      diagnostics = list(
        solver_status = overall_comp_status,
        objective_value = total_obj,
        n_variables = NA_integer_,
        solve_time = elapsed,
        optimality_gap = NA_real_
      )
    ))
  }
  # -- End component decomposition --------------------------------------------

  # -- spopt backend -----------------------------------------------------------
  if (solver == "spopt") {
    return(.solve_partition_spopt(
      full_sparse_distances = full_sparse_distances,
      tracts = tracts,
      D_max = D_max,
      max_workload_per_zone = max_workload_per_zone,
      strategy = effective_strategy,
      verbose = verbose
    ))
  }
  # -- End spopt backend -------------------------------------------------------

  # -- uncap_then_split strategy ------------------------------------------------
  if (effective_strategy == "uncap_then_split" && capacitated) {
    return(.uncap_then_split(
      full_sparse_distances = full_sparse_distances,
      tracts = tracts,
      K0 = K0,
      K_max = K_max,
      D_max = D_max,
      max_workload_per_zone = max_workload_per_zone,
      solver = solver,
      max_time = max_time,
      rel_tol = rel_tol,
      verbose = verbose
    ))
  }

  # -- End uncap_then_split ----------------------------------------------------

  for (K in seq(K0, K_max)) {
    cli::cli_alert("Trying K = {K}...")

    result <- surveyzones_solve_fixed_K(
      full_sparse_distances = full_sparse_distances,
      tracts = tracts,
      K = K,
      D_max = D_max,
      max_workload_per_zone = max_workload_per_zone,
      solver = solver,
      max_time = max_time,
      rel_tol = rel_tol,
      verbose = verbose
    )

    if (.is_solution_status(result$diagnostics$solver_status)) {
      cli::cli_alert_success(
        "Feasible at K = {K} (objective = {round(result$diagnostics$objective_value, 2)})"
      )
      return(result)
    }

    cli::cli_alert_warning("Infeasible at K = {K}, incrementing...")
  }

  cli::cli_alert_danger("No feasible solution found up to K_max = {K_max}.")
  # Return an infeasible result
  .empty_zone_solution(
    solver_status = "infeasible",
    n_variables = NA_integer_,
    solve_time = NA_real_
  )
}


#' Add Assignment Constraints to MILP
#'
#' Each tract must be assigned to exactly one center.
#'
#' @param model An ompr MILP model.
#' @param ps_by_i List mapping tract index to sparse-pair indices.
#' @param n Integer. Number of tracts.
#'
#' @return Updated MILP model with assignment constraints.
#' @keywords internal
.add_assignment_constraints <- function(model, ps_by_i, n) {
  for (i_val in seq_len(n)) {
    ps <- ps_by_i[[as.character(i_val)]]
    model <- model |>
      ompr::add_constraint(ompr::sum_over(x[p], p = ps) == 1)
  }
  model
}


#' Add Capacity Constraints to MILP
#'
#' Each center cannot exceed the maximum workload.
#'
#' @param model An ompr MILP model.
#' @param ps_by_j List mapping center index to sparse-pair indices.
#' @param pair_i Integer vector mapping sparse pair index to tract index.
#' @param w Numeric vector of tract workloads.
#' @param m Integer. Number of candidate centers.
#' @param max_workload_per_zone Numeric. Maximum workload per zone.
#'
#' @return Updated MILP model with capacity constraints.
#' @keywords internal
.add_capacity_constraints <- function(model, ps_by_j, pair_i, w, m, max_workload_per_zone) {
  for (j_val in seq_len(m)) {
    ps <- ps_by_j[[as.character(j_val)]]
    if (length(ps) > 0) {
      model <- model |>
        ompr::add_constraint(
          ompr::sum_over(w[pair_i[p]] * x[p], p = ps) <= max_workload_per_zone
        )
    }
  }
  model
}


#' Build Solver Control Parameters
#'
#' Constructs a list of solver-specific parameters for ompr.roi::with_ROI().
#' Different solvers accept different parameter names.
#'
#' @param solver Character scalar: `"glpk"`, `"highs"`, `"cbc"`, or `"symphony"`.
#' @param max_time Numeric. Maximum solve time in seconds.
#' @param rel_tol Numeric. Relative tolerance (0-1).
#' @param verbose Logical. Whether to show solver output.
#'
#' @return List of named parameters for with_ROI().
#' @keywords internal
.solver_control_params <- function(solver, max_time, rel_tol, verbose) {
  if (solver == "symphony") {
    list(
      solver = solver,
      max_time = as.numeric(max_time),
      gap_limit = rel_tol * 100,
      verbose = verbose
    )
  } else {
    list(
      solver = solver,
      max_time = as.numeric(max_time),
      rel_tol = rel_tol,
      verbose = verbose
    )
  }
}


#' Solve the Capacitated p-Median Problem for Fixed K
#'
#' Formulates and solves a mixed-integer linear program:
#' assign each tract to exactly one center, open exactly K centers,
#' respect workload capacity, and minimise total weighted distance.
#'
#' @inheritParams surveyzones_build_zones
#' @param full_sparse_distances Tibble of all distance pairs (unfiltered).
#'   Will be filtered internally to D_max.
#' @param K Integer.  Number of zones (centers) to open.
#'
#' @return A list with `assignments`, `zones`, and `diagnostics`.
#' @keywords internal
surveyzones_solve_fixed_K <- function(
  full_sparse_distances,
  tracts,
  K,
  D_max,
  max_workload_per_zone,
  solver = "glpk",
  max_time = 300,
  rel_tol = 0.01,
  verbose = FALSE
) {
  max_variables <- 500000L
  start_time <- proc.time()["elapsed"]

  tract_ids <- as.character(tracts$tract_id)
  n <- length(tract_ids)

  cli::cli_alert_info("  K={K}: {n} tracts, building MILP...")

  # Build workload lookup
  workload <- stats::setNames(tracts$expected_service_time, tract_ids)

  # Filter distances to this partition and D_max
  part_ids <- tract_ids
  sparse_distances <- full_sparse_distances |>
    dplyr::filter(
      distance <= D_max,
      origin_id %in% part_ids,
      destination_id %in% part_ids
    )

  # Candidate centers are all tracts
  m <- n

  # Sparse pairs: (tract index i, candidate index j) where distance exists
  # Map IDs to indices
  tract_to_idx <- stats::setNames(seq_len(n), tract_ids)
  cand_to_jdx <- stats::setNames(seq_len(m), tract_ids)

  dt <- sparse_distances

  # Also allow self-assignment (center assigned to itself with distance 0)
  self_pairs <- tibble::tibble(
    origin_id = tract_ids,
    destination_id = tract_ids,
    distance = 0
  )
  dt <- dplyr::bind_rows(dt, self_pairs) |>
    dplyr::distinct(origin_id, destination_id, .keep_all = TRUE)

  # Map to numeric indices
  dt <- dt |>
    dplyr::mutate(
      i = tract_to_idx[origin_id],
      j = cand_to_jdx[destination_id]
    ) |>
    dplyr::filter(!is.na(i), !is.na(j))

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
    "  K={K}: {n_x} x-variables, {m} potential centers, {n} assignment constraints, {n_cap} capacity constraints"
  )

  # Check that every tract can reach at least one candidate
  reachable <- unique(dt$i)
  unreachable <- setdiff(seq_len(n), reachable)
  if (length(unreachable) > 0) {
    unreach_ids <- tract_ids[unreachable]
    cli::cli_alert_warning(
      "{length(unreachable)} tract{?s} unreachable: {.val {head(unreach_ids, 5)}}{if (length(unreach_ids) > 5) '...'}"
    )
    return(.empty_zone_solution(
      solver_status = "infeasible",
      n_variables = as.integer(n_x),
      solve_time = as.numeric(proc.time()["elapsed"] - start_time)
    ))
  }

  # Build the MILP with ompr
  # Sparse x: only for pairs that exist in dt
  # Index the sparse pairs as 1..n_x
  pair_i <- dt$i
  pair_j <- dt$j
  pair_d <- dt$distance
  ps_by_i <- split(seq_len(n_x), pair_i)
  ps_by_j <- split(seq_len(n_x), pair_j)

  # Workload for each tract
  w <- workload[tract_ids]

  model <- ompr::MILPModel() |>
    ompr::add_variable(x[p], p = 1:n_x, type = "binary") |>
    ompr::add_variable(y[j], j = 1:m, type = "binary") |>
    ompr::set_objective(
      ompr::sum_over(pair_d[p] * x[p], p = 1:n_x),
      sense = "min"
    ) |>
    # x[p] <= y[j] -- can only assign to an open center
    ompr::add_constraint(
      x[p] <= y[pair_j[p]],
      p = 1:n_x
    ) |>
    # Exactly K centers
    ompr::add_constraint(
      ompr::sum_over(y[j], j = 1:m) == K
    )

  # Add assignment and capacity constraints
  model <- .add_assignment_constraints(model, ps_by_i, n)
  if (capacitated) {
    model <- .add_capacity_constraints(model, ps_by_j, pair_i, w, m, max_workload_per_zone)
  }

  # Solve with appropriate solver parameters
  solver_label <- toupper(solver)
  cli::cli_alert_info(
    "  K={K}: Solving MILP with {solver_label} (max_time={max_time}s, rel_tol={rel_tol})..."
  )
  solve_t0 <- proc.time()[["elapsed"]]
  solved <- .solve_model_with_status(
    model = model,
    solver = solver,
    max_time = max_time,
    rel_tol = rel_tol,
    verbose = verbose
  )
  result <- solved$result
  solve_elapsed <- proc.time()[["elapsed"]] - solve_t0

  solve_time <- as.numeric(proc.time()["elapsed"] - start_time)
  raw_status <- solved$raw_status
  status <- solved$status
  obj_try <- solved$objective_value
  if (status == "feasible" && raw_status != "feasible") {
    cli::cli_alert_warning(
      "  K={K}: solver returned {raw_status} but produced a feasible incumbent (obj={round(obj_try, 2)})"
    )
  }

  obj_val <- tryCatch(round(ompr::objective_value(result), 2), error = function(e) NA_real_)
  cli::cli_alert_info(
    "  K={K}: solver returned {.val {status}} in {round(solve_elapsed, 1)}s, objective={obj_val}"
  )

  if (!.is_solution_status(status)) {
    return(.empty_zone_solution(
      solver_status = status,
      n_variables = as.integer(n_x),
      solve_time = solve_time
    ))
  }

  # Extract solution
  x_sol <- ompr::get_solution(result, x[p])
  active <- x_sol[x_sol$value > 0.5, ]

  cli::cli_alert_info(
    "  K={K}: {nrow(active)} active assignments, objective = {round(ompr::objective_value(result), 2)}"
  )

  # Validate all tracts are assigned
  if (nrow(active) < n) {
    cli::cli_alert_warning(
      "  K={K}: only {nrow(active)}/{n} tracts assigned, treating as infeasible"
    )
    return(.empty_zone_solution(
      solver_status = "infeasible",
      n_variables = as.integer(n_x),
      solve_time = solve_time
    ))
  }

  # Build assignments with all columns in dplyr pipeline
  assignments <- tibble::tibble(
    tract_id = tract_ids[pair_i[active$p]],
    center_id = tract_ids[pair_j[active$p]],
    distance_to_center = pair_d[active$p]
  ) |>
    dplyr::mutate(
      zone_id = center_id,
      partition_id = NA_character_
    )

  # Zone summary with all columns in one pipeline
  zones <- assignments |>
    dplyr::summarise(
      center_tract_id = unique(center_id),
      total_workload = sum(workload[tract_id]),
      n_tracts = dplyr::n(),
      partition_id = NA_character_,
      .by = zone_id
    )

  # Compute zone diameters using helper function
  zones$diameter <- purrr::map_dbl(
    zones$zone_id,
    \(zid) .compute_zone_diameter(
      assignments$tract_id[assignments$zone_id == zid],
      full_sparse_distances
    )
  )

  cli::cli_alert_info(
    "  K={K}: {nrow(zones)} zones, workload range [{round(min(zones$total_workload), 1)}-{round(max(zones$total_workload), 1)}], max diameter {round(max(zones$diameter, na.rm = TRUE), 2)} (distance units)"
  )

  list(
    assignments = assignments |>
      dplyr::mutate(group_id = zone_id) |>
      dplyr::select(tract_id, zone_id, partition_id, center_id, distance_to_center, group_id),
    zones = zones |>
      dplyr::mutate(group_id = zone_id) |>
      dplyr::select(zone_id, partition_id, center_tract_id, total_workload, diameter, n_tracts, group_id),
    diagnostics = list(
      solver_status = status,
      objective_value = ompr::objective_value(result),
      n_variables = as.integer(n_x),
      solve_time = solve_time,
      optimality_gap = solved$optimality_gap
    )
  )
}


#' Uncapacitated-then-Split Heuristic
#'
#' Solves the uncapacitated p-median (fast), then recursively splits any zone
#' whose total workload exceeds `max_workload_per_zone` into sub-zones.
#' Each split sub-problem contains only the tracts of the oversized zone, so
#' it is tiny and solved almost instantly.
#'
#' @inheritParams surveyzones_build_zones
#' @param K0 Integer. Initial number of zones (lower bound).
#' @param K_max Integer. Upper bound on number of zones.
#' @return A list with `assignments`, `zones`, and `diagnostics`.
#' @keywords internal
.uncap_then_split <- function(
  full_sparse_distances,
  tracts,
  K0,
  K_max,
  D_max,
  max_workload_per_zone,
  solver,
  max_time,
  rel_tol,
  verbose
) {
  start_time <- proc.time()[["elapsed"]]

  cli::cli_alert_info(
    "Strategy: uncap_then_split \u2014 solving uncapacitated from K = {K0} to K_max = {K_max}, then splitting oversized zones"
  )

  # Phase 1: uncapacitated solve — increment K until feasible
  uncap_result <- NULL
  for (K in seq(K0, K_max)) {
    cli::cli_alert("Trying uncapacitated K = {K}...")
    uncap_result <- surveyzones_solve_fixed_K(
      full_sparse_distances = full_sparse_distances,
      tracts = tracts,
      K = K,
      D_max = D_max,
      max_workload_per_zone = Inf,
      solver = solver,
      max_time = max_time,
      rel_tol = rel_tol,
      verbose = verbose
    )
    if (.is_solution_status(uncap_result$diagnostics$solver_status)) {
      cli::cli_alert_success("Uncapacitated solve feasible at K = {K}")
      break
    }
    cli::cli_alert_warning("Uncapacitated K = {K} failed (status: {uncap_result$diagnostics$solver_status})")
  }

  if (!.is_solution_status(uncap_result$diagnostics$solver_status)) {
    cli::cli_alert_danger("Uncapacitated solve failed for all K in {K0}..{K_max}")
    return(uncap_result)
  }
  overall_status <- uncap_result$diagnostics$solver_status

  # Phase 2: find oversized zones and split them
  workload_by_tract <- stats::setNames(
    tracts$expected_service_time, as.character(tracts$tract_id)
  )

  zone_summary <- uncap_result$assignments |>
    dplyr::mutate(.wl = workload_by_tract[.data$tract_id]) |>
    dplyr::summarise(total_wl = sum(.data$.wl), .by = "zone_id")

  oversized <- zone_summary |>
    dplyr::filter(.data$total_wl > max_workload_per_zone) |>
    dplyr::pull(.data$zone_id)

  n_over <- length(oversized)
  if (n_over == 0L) {
    cli::cli_alert_success("All zones satisfy the workload cap \u2014 no splitting needed")
    # Each zone is its own group
    uncap_result$assignments$group_id <- uncap_result$assignments$zone_id
    uncap_result$zones$group_id <- uncap_result$zones$zone_id
    return(uncap_result)
  }

  cli::cli_alert_info("{n_over} oversized zone{?s} to split")

  good_assignments <- uncap_result$assignments |>
    dplyr::filter(!.data$zone_id %in% oversized) |>
    dplyr::mutate(group_id = .data$zone_id)
  good_zones <- uncap_result$zones |>
    dplyr::filter(!.data$zone_id %in% oversized) |>
    dplyr::mutate(group_id = .data$zone_id)

  split_assignments <- vector("list", n_over)
  split_zones <- vector("list", n_over)

  for (idx in seq_along(oversized)) {
    zid <- oversized[[idx]]
    zone_tract_ids <- uncap_result$assignments$tract_id[
      uncap_result$assignments$zone_id == zid
    ]
    zone_tracts <- tracts |> dplyr::filter(.data$tract_id %in% zone_tract_ids)
    zone_wl <- sum(workload_by_tract[zone_tract_ids])
    K_split <- ceiling(zone_wl / max_workload_per_zone)

    n_tracts_in_zone <- length(zone_tract_ids)
    if (K_split > n_tracts_in_zone) {
      cli::cli_alert_danger(
        "  Zone {zid} cannot be split under cap: minimum K={K_split} exceeds tract count {n_tracts_in_zone}"
      )
      return(.empty_zone_solution(
        solver_status = "infeasible",
        n_variables = uncap_result$diagnostics$n_variables,
        solve_time = proc.time()[["elapsed"]] - start_time
      ))
    }

    cli::cli_alert_info(
      "  Splitting zone {zid} ({n_tracts_in_zone} tracts, wl={round(zone_wl,1)}) into K={K_split}"
    )

    sub_dist <- full_sparse_distances |>
      dplyr::filter(
        .data$origin_id %in% zone_tract_ids,
        .data$destination_id %in% zone_tract_ids
      )

    # Try increasing K until feasible (at K = n_tracts, trivially feasible)
    sub_result <- NULL
    for (K_try in seq.int(K_split, n_tracts_in_zone)) {
      if (K_try > K_split) {
        cli::cli_alert_info(
          "  Retrying zone {zid} with K={K_try}"
        )
      }
      sub_result <- surveyzones_solve_fixed_K(
        full_sparse_distances = sub_dist,
        tracts = zone_tracts,
        K = K_try,
        D_max = D_max,
        max_workload_per_zone = max_workload_per_zone,
        solver = solver,
        max_time = max_time,
        rel_tol = rel_tol,
        verbose = verbose
      )
      if (.is_solution_status(sub_result$diagnostics$solver_status)) break
    }

    if (!.is_solution_status(sub_result$diagnostics$solver_status)) {
      cli::cli_alert_danger(
        "  Split of zone {zid} failed after trying K={K_split}..{n_tracts_in_zone}; cannot satisfy workload cap"
      )
      return(.empty_zone_solution(
        solver_status = "infeasible",
        n_variables = uncap_result$diagnostics$n_variables,
        solve_time = proc.time()[["elapsed"]] - start_time
      ))
    } else {
      if (sub_result$diagnostics$solver_status == "feasible") {
        overall_status <- "feasible"
      }
      if (K_try > K_split) {
        cli::cli_alert_success(
          "  Zone {zid} split succeeded at K={K_try} (minimum K={K_split} was infeasible)"
        )
      }
      # Tag split sub-zones with parent phase-1 zone as group_id
      split_assignments[[idx]] <- sub_result$assignments |>
        dplyr::mutate(group_id = zid)
      split_zones[[idx]] <- sub_result$zones |>
        dplyr::mutate(group_id = zid)
    }
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
    "uncap_then_split: {nrow(final_zones)} zones in {round(elapsed, 1)}s"
  )

  list(
    assignments = final_assignments,
    zones = final_zones,
    diagnostics = list(
      solver_status = overall_status,
      objective_value = sum(final_assignments$distance_to_center, na.rm = TRUE),
      n_variables = uncap_result$diagnostics$n_variables,
      solve_time = elapsed,
      optimality_gap = NA_real_
    )
  )
}


#' @keywords internal
.symmetrize_distances_max <- function(distance_table) {
  if (nrow(distance_table) == 0L) {
    return(tibble::tibble(
      origin_id = character(0),
      destination_id = character(0),
      distance = numeric(0)
    ))
  }

  undirected <- distance_table |>
    dplyr::transmute(
      origin_id = as.character(.data$origin_id),
      destination_id = as.character(.data$destination_id),
      distance = .data$distance
    ) |>
    dplyr::filter(.data$origin_id != .data$destination_id) |>
    dplyr::mutate(
      a = pmin(.data$origin_id, .data$destination_id),
      b = pmax(.data$origin_id, .data$destination_id)
    ) |>
    dplyr::summarise(
      distance = {
        d <- .data$distance[is.finite(.data$distance)]
        if (length(d) == 0L) NA_real_ else max(d)
      },
      .by = c("a", "b")
    ) |>
    dplyr::filter(is.finite(.data$distance))

  dplyr::bind_rows(
    undirected |>
      dplyr::transmute(origin_id = .data$a, destination_id = .data$b, distance = .data$distance),
    undirected |>
      dplyr::transmute(origin_id = .data$b, destination_id = .data$a, distance = .data$distance)
  ) |>
    dplyr::distinct(.data$origin_id, .data$destination_id, .keep_all = TRUE)
}






#' Compute Zone Diameter
#'
#' Calculates the maximum pairwise distance within a zone.
#' Uses the full (unfiltered) sparse distances so peripheral pairs
#' beyond D_max are included. Returns NA when any member pair is
#' missing from the distance table (true diameter unknown).
#'
#' @param tract_ids Character vector of tracts in the zone.
#' @param full_sparse_distances Tibble of all distance pairs.
#'
#' @return Numeric scalar: maximum distance, or NA if incomplete.
#' @keywords internal
.compute_zone_diameter <- function(tract_ids, full_sparse_distances) {
  n_members <- length(tract_ids)
  if (n_members <= 1L) return(0)

  intra <- full_sparse_distances |>
    dplyr::filter(
      origin_id %in% tract_ids,
      destination_id %in% tract_ids,
      origin_id != destination_id
    )

  if (nrow(intra) == 0L) return(NA_real_)

  # Check completeness: need all n*(n-1)/2 unique unordered pairs
  pairs <- unique(paste0(
    pmin(intra$origin_id, intra$destination_id), "|",
    pmax(intra$origin_id, intra$destination_id)
  ))
  expected <- n_members * (n_members - 1L) %/% 2L
  if (length(pairs) < expected) return(NA_real_)

  max(intra$distance)
}


#' Validate and Load a Solver Plugin
#'
#' Checks that the requested solver is supported and that the
#' corresponding ROI plugin package (or spopt) is available.
#'
#' @param solver Character scalar: `"glpk"`, `"highs"`, `"cbc"`,
#'   `"symphony"`, or `"spopt"`.
#' @param call Caller environment for error reporting.
#' @return `solver`, invisibly.
#' @keywords internal
validate_solver <- function(solver, call = rlang::caller_env()) {
  supported <- c("glpk", "highs", "cbc", "symphony", "spopt")
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
        "i" = "For MILP solvers, install the matching ROI plugin: {.pkg ROI.plugin.{solver}}"
      ),
      call = call
    )
  }
  if (solver == "spopt") {
    rlang::check_installed("spopt", reason = "to use the spopt solver backend.")
    return(invisible(solver))
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

  if (anyNA(tracts$expected_service_time)) {
    cli::cli_abort(
      "{.field expected_service_time} must not contain NA values.",
      call = call
    )
  }

  if (any(!is.finite(tracts$expected_service_time))) {
    cli::cli_abort(
      "{.field expected_service_time} must be finite for all tracts.",
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

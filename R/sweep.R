#' Parameter Sweep Across D_max and Workload Combinations
#'
#' Runs [surveyzones_build_zones()] for each combination of `D_max`
#' and `max_workload_per_zone`, collecting summary diagnostics.
#'
#' @param sparse_distances Sparse distance table.
#' @param tracts Tract table.
#' @param D_max_values Numeric vector of D_max values to try.
#' @param max_workload_values Numeric vector of max_workload_per_zone
#'   values to try.
#' @param enforce_partition Logical.
#' @param solver Character scalar.  Which MILP solver backend to use.
#'   One of `"glpk"` (default), `"highs"`, or `"cbc"`.
#' @param ... Additional arguments passed to [surveyzones_build_zones()].
#'
#' @return A tibble with one row per parameter combination and columns:
#'   `D_max`, `max_workload_per_zone`, `K_total`,
#'   `n_infeasible_partitions`, `mean_diameter`, `max_diameter`,
#'   `mean_workload`, `objective_total`.
#'
#' @export
surveyzones_sweep <- function(
    sparse_distances,
    tracts,
    D_max_values,
    max_workload_values,
    enforce_partition = TRUE,
    solver = "glpk",
    ...) {

  grid <- expand.grid(
    D_max = D_max_values,
    max_workload_per_zone = max_workload_values,
    stringsAsFactors = FALSE
  )

  cli::cli_alert_info(
    "Running {nrow(grid)} parameter combination{?s}..."
  )

  results <- vector("list", nrow(grid))

  for (i in seq_len(nrow(grid))) {
    d_max <- grid$D_max[i]
    mw <- grid$max_workload_per_zone[i]

    cli::cli_h3("D_max = {d_max}, max_workload = {mw}")

    plan <- tryCatch(
      surveyzones_build_zones(
        sparse_distances = sparse_distances,
        tracts = tracts,
        D_max = d_max,
        max_workload_per_zone = mw,
        enforce_partition = enforce_partition,
        solver = solver,
        ...
      ),
      error = function(e) {
        cli::cli_alert_danger("Error: {conditionMessage(e)}")
        NULL
      }
    )

    if (is.null(plan) || nrow(plan$zones) == 0) {
      results[[i]] <- tibble::tibble(
        D_max = d_max,
        max_workload_per_zone = mw,
        K_total = NA_integer_,
        n_infeasible_partitions = NA_integer_,
        mean_diameter = NA_real_,
        max_diameter = NA_real_,
        mean_workload = NA_real_,
        objective_total = NA_real_
      )
    } else {
      n_infeasible <- sum(plan$diagnostics$solver_status != "optimal")
      results[[i]] <- tibble::tibble(
        D_max = d_max,
        max_workload_per_zone = mw,
        K_total = nrow(plan$zones),
        n_infeasible_partitions = n_infeasible,
        mean_diameter = mean(plan$zones$diameter, na.rm = TRUE),
        max_diameter = max(plan$zones$diameter, na.rm = TRUE),
        mean_workload = mean(plan$zones$total_workload),
        objective_total = sum(plan$diagnostics$objective_value, na.rm = TRUE)
      )
    }
  }

  do.call(rbind, results)
}


#' Compare Multiple Plans
#'
#' Extracts summary diagnostics from a named list of
#' `surveyzones_plan` objects and returns a comparison tibble.
#'
#' @param plans A named list of `surveyzones_plan` objects.
#'
#' @return A tibble with one row per plan and columns:
#'   `plan_name`, `n_partitions`, `n_tracts`, `n_zones`,
#'   `workload_min`, `workload_max`, `workload_mean`,
#'   `diameter_min`, `diameter_max`, `diameter_mean`,
#'   `objective_total`.
#'
#' @export
surveyzones_compare <- function(plans) {
  if (!is.list(plans) || is.null(names(plans))) {
    cli::cli_abort("{.arg plans} must be a named list of {.cls surveyzones_plan} objects.")
  }

  rows <- lapply(names(plans), function(nm) {
    plan <- plans[[nm]]
    if (!inherits(plan, "surveyzones_plan")) {
      cli::cli_abort("Element {.val {nm}} is not a {.cls surveyzones_plan}.")
    }
    z <- plan$zones
    tibble::tibble(
      plan_name = nm,
      n_partitions = length(unique(z$partition_id)),
      n_tracts = nrow(plan$assignments),
      n_zones = nrow(z),
      workload_min = if (nrow(z) > 0) min(z$total_workload) else NA_real_,
      workload_max = if (nrow(z) > 0) max(z$total_workload) else NA_real_,
      workload_mean = if (nrow(z) > 0) mean(z$total_workload) else NA_real_,
      diameter_min = if (nrow(z) > 0) min(z$diameter, na.rm = TRUE) else NA_real_,
      diameter_max = if (nrow(z) > 0) max(z$diameter, na.rm = TRUE) else NA_real_,
      diameter_mean = if (nrow(z) > 0) mean(z$diameter, na.rm = TRUE) else NA_real_,
      objective_total = sum(plan$diagnostics$objective_value, na.rm = TRUE)
    )
  })

  result <- do.call(rbind, rows)
  result
}

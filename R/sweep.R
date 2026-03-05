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

#' Zone Workload Summary
#'
#' Aggregate `expected_service_time` by zone.
#'
#' @param plan A `surveyzones_plan` object.
#' @param tracts The tract table used to build the plan.
#'
#' @return A tibble with columns `zone_id`, `partition_id`,
#'   `total_workload`, `n_tracts`.
#'
#' @export
surveyzones_zone_workload <- function(plan, tracts) {
  if (!inherits(plan, "surveyzones_plan")) {
    cli::cli_abort("{.arg plan} must be a {.cls surveyzones_plan} object.")
  }

  workload <- stats::setNames(tracts$expected_service_time,
                               as.character(tracts$tract_id))

  plan$assignments |>
    dplyr::summarise(
      total_workload = sum(workload[tract_id]),
      n_tracts = dplyr::n(),
      .by = c(zone_id, partition_id)
    )
}

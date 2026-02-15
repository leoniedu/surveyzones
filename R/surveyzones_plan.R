#' Create a surveyzones_plan Object
#'
#' Constructor for the S3 class returned by [surveyzones_build_zones()].
#'
#' @param assignments A tibble with columns `tract_id`, `zone_id`,
#'   `partition_id`, `center_id`, `travel_time_to_center`.
#' @param zones A tibble with columns `zone_id`, `partition_id`,
#'   `center_tract_id`, `total_workload`, `diameter`, `n_tracts`.
#' @param parameters A named list of model parameters.
#' @param diagnostics A named list of solver diagnostics.
#'
#' @return A `surveyzones_plan` object.
#' @keywords internal
new_surveyzones_plan <- function(assignments, zones, parameters,
                                 diagnostics) {
  structure(
    list(
      assignments = assignments,
      zones = zones,
      parameters = parameters,
      diagnostics = diagnostics,
      sequence = NULL
    ),
    class = "surveyzones_plan"
  )
}


#' @export
print.surveyzones_plan <- function(x, ...) {
  n_partitions <- length(unique(x$zones$partition_id))
  n_tracts <- nrow(x$assignments)
  n_zones <- nrow(x$zones)

  cli::cli_h1("surveyzones plan")
  cli::cli_bullets(c(
    "*" = "{n_partitions} partition{?s}, {n_tracts} tract{?s}, {n_zones} zone{?s}",
    "*" = "D_max = {x$parameters$D_max}, max_workload = {if (is.infinite(x$parameters$max_workload_per_zone)) 'Inf' else x$parameters$max_workload_per_zone}{if (!is.null(x$parameters$target_zone_size)) paste0(', target_zone_size = ', x$parameters$target_zone_size) else ''}"
  ))

  if (n_zones > 0) {
    wl <- x$zones$total_workload
    dm <- x$zones$diameter
    cli::cli_bullets(c(
      "*" = "Workload per zone: {round(min(wl), 1)}\u2013{round(max(wl), 1)} (mean {round(mean(wl), 1)})",
      "*" = "Diameter per zone: {round(min(dm, na.rm = TRUE), 2)}\u2013{round(max(dm, na.rm = TRUE), 2)} (mean {round(mean(dm, na.rm = TRUE), 2)})"
    ))
  }

  # Solver status
  statuses <- unique(x$diagnostics$solver_status)
  if (all(statuses == "optimal")) {
    cli::cli_bullets(c("v" = "All partitions solved optimally"))
  } else {
    n_fail <- sum(statuses != "optimal")
    cli::cli_bullets(c("!" = "{n_fail} partition{?s} not optimal: {.val {statuses}}"))
  }

  has_seq <- !is.null(x$sequence)
  if (has_seq) {
    cli::cli_bullets(c("*" = "Sequencing: computed"))
  }

  invisible(x)
}


#' @export
summary.surveyzones_plan <- function(object, ...) {
  zones <- object$zones

  if (nrow(zones) == 0) {
    cli::cli_alert_warning("No zones in this plan.")
    return(invisible(object))
  }

  # Per-partition summary
  summary_tbl <- zones |>
    tidyr::nest(data = -partition_id) |>
    dplyr::mutate(
      K = purrr::map_int(data, nrow),
      n_tracts = purrr::map_dbl(data, \(z) sum(z$n_tracts)),
      workload_min = purrr::map_dbl(data, \(z) min(z$total_workload)),
      workload_max = purrr::map_dbl(data, \(z) max(z$total_workload)),
      workload_mean = purrr::map_dbl(data, \(z) mean(z$total_workload)),
      diameter_min = purrr::map_dbl(data, \(z) min(z$diameter, na.rm = TRUE)),
      diameter_max = purrr::map_dbl(data, \(z) max(z$diameter, na.rm = TRUE)),
      diameter_mean = purrr::map_dbl(data, \(z) mean(z$diameter, na.rm = TRUE))
    ) |>
    dplyr::select(-data)

  cli::cli_h1("surveyzones plan summary")
  print(summary_tbl)

  # Overall diagnostics
  total_obj <- sum(object$diagnostics$objective_value, na.rm = TRUE)
  total_time <- sum(object$diagnostics$solve_time, na.rm = TRUE)
  cli::cli_bullets(c(
    "*" = "Total objective: {round(total_obj, 2)}",
    "*" = "Total solve time: {round(total_time, 1)}s"
  ))

  invisible(summary_tbl)
}

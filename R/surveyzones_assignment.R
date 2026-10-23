#' Create a surveyzones_assignment Object
#'
#' Constructor for the S3 class returned by [surveyzones_assign_teams()].
#'
#' @param assignments Tibble with `tract_id`, `team_id`, `roteiro_id`, `mode`, `cost`.
#' @param roteiros Tibble with `roteiro_id`, `team_id`, `center_point_id`,
#'   `center_tract_id`, `n_tracts`, `workload`, `span_m`, `cost`.
#' @param loads Tibble with `team_id`, `workload`, `n_roteiros`, `n_loose`,
#'   `loose_cost`, `lower`, `upper`.
#' @param origins Tibble with `team_id`, `origin_tract_id`.
#' @param parameters Named list of model parameters.
#' @param diagnostics Named list of solver diagnostics.
#'
#' @return A `surveyzones_assignment` object.
#' @keywords internal
new_surveyzones_assignment <- function(assignments, roteiros, loads, origins,
                                       parameters, diagnostics) {
  structure(
    list(
      assignments = assignments,
      roteiros = roteiros,
      loads = loads,
      origins = origins,
      parameters = parameters,
      diagnostics = diagnostics
    ),
    class = "surveyzones_assignment"
  )
}

#' @export
print.surveyzones_assignment <- function(x, ...) {
  cli::cli_h1("surveyzones assignment")
  status <- x$diagnostics$solver_status
  if (!is.null(x$diagnostics$elastic)) {
    cli::cli_bullets(c("x" = "INFEASIBLE \u2014 see $diagnostics$elastic"))
  }
  if (nrow(x$assignments) == 0L) {
    cli::cli_bullets(c("!" = "No assignment (solver status: {.val {status}})"))
    return(invisible(x))
  }
  modes <- table(x$assignments$mode)
  k_txt <- if (isTRUE(x$parameters$k_auto)) {
    "{x$parameters$k} (auto; minimum {x$parameters$k_min})"
  } else {
    "{x$parameters$k}"
  }
  if (!is.null(x$parameters$k)) {
    cli::cli_bullets(c("*" = paste0("k = ", k_txt)))
  }
  cli::cli_bullets(c(
    "*" = "{nrow(x$assignments)} tract{?s} to {nrow(x$loads)} team{?s}: {sum(modes[c('roteiro')], na.rm = TRUE)} in roteiros, {sum(modes[c('loose')], na.rm = TRUE)} loose, {sum(modes[c('fixed')], na.rm = TRUE)} fixed",
    "*" = "Workload per team: {round(min(x$loads$workload), 1)} to {round(max(x$loads$workload), 1)} (band {round(x$loads$lower[1], 1)} to {round(x$loads$upper[1], 1)})",
    "*" = "Objective: {round(x$diagnostics$objective_value, 2)}"
  ))
  if (status == "optimal") {
    cli::cli_bullets(c("v" = "Solved optimally"))
  } else {
    gap <- x$diagnostics$optimality_gap
    if (identical(status, "feasible") && length(gap) == 1L && is.finite(gap)) {
      cli::cli_bullets(c("!" = "Solver status: {.val {status}} (optimality gap {round(100 * gap, 2)}%)"))
    } else {
      cli::cli_bullets(c("!" = "Solver status: {.val {status}}"))
    }
  }
  invisible(x)
}

#' @export
summary.surveyzones_assignment <- function(object, ...) {
  if (nrow(object$assignments) == 0L) {
    cli::cli_alert_warning("No assignment in this object.")
    return(invisible(object$loads))
  }
  cli::cli_h1("surveyzones assignment summary")
  print(object$loads)
  cli::cli_bullets(c(
    "*" = "Total objective: {round(object$diagnostics$objective_value, 2)}",
    "*" = "Solve time: {round(object$diagnostics$solve_time, 1)}s"
  ))
  invisible(object$loads)
}

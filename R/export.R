#' Export a Zone Plan to Parquet
#'
#' Writes the complete zone plan to a directory of Parquet files
#' for archival and auditing.  Requires the arrow package.
#'
#' @param plan A `surveyzones_plan` object.
#' @param path Directory path where Parquet files will be written.
#'
#' @return `path`, invisibly.
#'
#' @export
surveyzones_export_plan <- function(plan, path) {
  if (!inherits(plan, "surveyzones_plan")) {
    cli::cli_abort("{.arg plan} must be a {.cls surveyzones_plan} object.")
  }

  rlang::check_installed("arrow", reason = "to export plans to Parquet")

  if (!dir.exists(path)) {
    dir.create(path, recursive = TRUE)
  }

  # Assignments
  arrow::write_parquet(
    plan$assignments,
    file.path(path, "assignments.parquet")
  )

  # Zones
  arrow::write_parquet(
    plan$zones,
    file.path(path, "zones.parquet")
  )

  # Sequence (if available)
  if (!is.null(plan$sequence)) {
    arrow::write_parquet(
      plan$sequence,
      file.path(path, "sequence.parquet")
    )
  }

  # Parameters as JSON-like tibble
  params_tbl <- tibble::tibble(
    parameter = names(plan$parameters),
    value = vapply(plan$parameters, function(v) {
      if (inherits(v, "POSIXct")) format(v) else as.character(v)
    }, character(1))
  )
  arrow::write_parquet(params_tbl, file.path(path, "parameters.parquet"))

  # Diagnostics
  diag_tbl <- tibble::as_tibble(plan$diagnostics)
  arrow::write_parquet(diag_tbl, file.path(path, "diagnostics.parquet"))

  cli::cli_alert_success("Plan exported to {.path {path}}.")
  invisible(path)
}

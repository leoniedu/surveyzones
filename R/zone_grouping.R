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

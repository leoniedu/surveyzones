#' Partition Tracts by Jurisdiction
#'
#' Splits a tract table into sub-problems by `partition_id`.
#' When `enforce_partition = TRUE` (the default) and a `partition_id`
#' column exists, the solver runs independently on each partition.
#'
#' @param tracts A data.frame or tibble with tract data.
#'   Must contain `tract_id`.  May contain `partition_id`.
#' @param enforce_partition Logical.  If `TRUE` (default), split
#'   by `partition_id`.  If `FALSE`, ignore partitions and return
#'   all tracts as a single group.
#'
#' @return A named list of data.frames, one per partition.
#'   If not partitioned, a single-element list named `"all"`.
#'
#' @export
surveyzones_partition <- function(tracts, enforce_partition = TRUE) {
  if (!is.data.frame(tracts)) {
    cli::cli_abort("{.arg tracts} must be a data.frame or tibble.")
  }

  if (!"tract_id" %in% names(tracts)) {
    cli::cli_abort("{.arg tracts} must contain a {.field tract_id} column.")
  }

  has_partition <- "partition_id" %in% names(tracts)

  if (enforce_partition && has_partition) {
    parts <- split(tracts, tracts$partition_id)
    cli::cli_alert_info(
      "Partitioned into {length(parts)} group{?s} by {.field partition_id}."
    )
    parts
  } else {
    if (enforce_partition && !has_partition) {
      cli::cli_alert_info(
        "No {.field partition_id} column found; solving as a single group."
      )
    }
    list(all = tracts)
  }
}

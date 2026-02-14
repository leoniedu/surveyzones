#' Write Sparse Distance Store
#'
#' Persists a sparse distance table to disk.
#'
#' @param distances A `data.table` as returned by
#'   [surveyzones_compute_sparse_distances()].
#' @param path File path (directory for arrow, file for others).
#' @param backend One of `"data.table"` (default), `"arrow"`, or
#'   `"duckdb"`.
#'
#' @return `path`, invisibly.
#'
#' @export
surveyzones_write_distance_store <- function(
    distances,
    path,
    backend = c("data.table", "arrow", "duckdb")) {

  backend <- match.arg(backend)

  switch(backend,
    data.table = {
      data.table::fwrite(distances, path)
    },
    arrow = {
      rlang::check_installed("arrow", reason = "to use the arrow backend")
      arrow::write_parquet(distances, path)
    },
    duckdb = {
      rlang::check_installed(c("duckdb", "DBI"), reason = "to use the duckdb backend")
      con <- DBI::dbConnect(duckdb::duckdb(), dbdir = path)
      on.exit(DBI::dbDisconnect(con, shutdown = TRUE), add = TRUE)
      DBI::dbWriteTable(con, "distances", distances, overwrite = TRUE)
    }
  )

  cli::cli_alert_success("Distances written to {.path {path}} ({backend}).")
  invisible(path)
}


#' Read Sparse Distance Store
#'
#' Reads a sparse distance table from disk.
#'
#' @param path File path.
#' @param backend One of `"data.table"` (default), `"arrow"`, or
#'   `"duckdb"`.
#'
#' @return A `data.table` with columns `origin_id`, `destination_id`,
#'   `travel_time`.
#'
#' @export
surveyzones_read_distance_store <- function(
    path,
    backend = c("data.table", "arrow", "duckdb")) {

  backend <- match.arg(backend)

  dt <- switch(backend,
    data.table = {
      data.table::fread(path)
    },
    arrow = {
      rlang::check_installed("arrow", reason = "to use the arrow backend")
      data.table::as.data.table(arrow::read_parquet(path))
    },
    duckdb = {
      rlang::check_installed(c("duckdb", "DBI"), reason = "to use the duckdb backend")
      con <- DBI::dbConnect(duckdb::duckdb(), dbdir = path, read_only = TRUE)
      on.exit(DBI::dbDisconnect(con, shutdown = TRUE), add = TRUE)
      data.table::as.data.table(DBI::dbReadTable(con, "distances"))
    }
  )

  data.table::setkey(dt, origin_id, destination_id)
  dt
}

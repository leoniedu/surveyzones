#' Clear the surveyzones disk cache
#'
#' Removes all cached results stored by [surveyzones_build_zones()].
#' Use this when inputs have changed structurally (e.g. different distance
#' engine or tract data) and you want to force a clean recompute.
#'
#' @param force Logical.  If `FALSE` (default), prints the cache size and asks
#'   for confirmation before deleting.  Set `TRUE` to skip the prompt.
#'
#' @return `TRUE` invisibly if the cache was cleared, `FALSE` if aborted.
#'
#' @export
surveyzones_clear_cache <- function(force = FALSE) {
  cache_dir <- .pkg$cache_dir
  if (!dir.exists(cache_dir)) {
    cli::cli_alert_info("Cache directory not found — nothing to clear.")
    return(invisible(FALSE))
  }

  files <- list.files(cache_dir, recursive = TRUE, full.names = TRUE)
  size  <- format(
    structure(sum(file.size(files)), class = "object_size"),
    units = "auto"
  )

  if (!force) {
    answer <- readline(
      sprintf("Clear surveyzones cache (%s)? [y/N] ", size)
    )
    if (!tolower(answer) %in% c("y", "yes")) {
      cli::cli_alert_warning("Aborted — cache not cleared.")
      return(invisible(FALSE))
    }
  }

  memoise::forget(surveyzones_build_zones_mem)
  cli::cli_alert_success("Cache cleared ({size} freed).")
  invisible(TRUE)
}

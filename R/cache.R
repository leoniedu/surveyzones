#' Clear the surveyzones disk cache
#'
#' Removes all cached results stored by [surveyzones_build_zones()],
#' [surveyzones_assign_teams()] and the Google Routes API engine ([surveyzones_engine_google()]).
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
  cache_dirs <- c(.pkg$cache_dir, .pkg$assign_cache_dir)
  cache_dirs <- cache_dirs[dir.exists(cache_dirs)]
  if (length(cache_dirs) == 0L) {
    cli::cli_alert_info("Cache directory not found \u2014 nothing to clear.")
    return(invisible(FALSE))
  }

  files <- list.files(cache_dirs, recursive = TRUE, full.names = TRUE)
  size  <- format(
    structure(sum(file.size(files)), class = "object_size"),
    units = "auto"
  )

  if (!force) {
    answer <- readline(
      sprintf("Clear surveyzones cache (%s)? [y/N] ", size)
    )
    if (!tolower(answer) %in% c("y", "yes")) {
      cli::cli_alert_warning("Aborted \u2014 cache not cleared.")
      return(invisible(FALSE))
    }
  }

  memoise::forget(surveyzones_build_zones_single_mem)
  memoise::forget(.assign_teams_core_mem)
  # keep the directories: cachem treats a missing directory as a destroyed cache
  unlink(list.files(cache_dirs, full.names = TRUE, all.files = TRUE, no.. = TRUE), recursive = TRUE)
  cli::cli_alert_success("Cache cleared ({size} freed).")
  invisible(TRUE)
}

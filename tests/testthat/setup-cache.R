# Redirect the disk-memoised solvers to a temp dir for the whole test run so
# tests never read or write (or clear) the user's cache. .onLoad assigns them
# with `<<-`, which lands in globalenv.
local({
  tmp <- tempfile("surveyzones-test-cache-")
  dir.create(tmp)
  pkg <- surveyzones:::.pkg
  old <- list(
    build = get("surveyzones_build_zones_single_mem", envir = globalenv()),
    assign = get(".assign_teams_core_mem", envir = globalenv()),
    dirs = as.list(pkg)[c("cache_dir", "assign_cache_dir")]
  )
  assign("surveyzones_build_zones_single_mem", memoise::memoise(
    surveyzones:::surveyzones_build_zones_single,
    cache = cachem::cache_disk(dir = file.path(tmp, "build_zones"))
  ), envir = globalenv())
  assign(".assign_teams_core_mem", memoise::memoise(
    surveyzones:::.assign_teams_core,
    omit_args = "verbose",
    cache = cachem::cache_disk(dir = file.path(tmp, "assign_teams"))
  ), envir = globalenv())
  pkg$cache_dir <- file.path(tmp, "build_zones")
  pkg$assign_cache_dir <- file.path(tmp, "assign_teams")
  withr::defer({
    assign("surveyzones_build_zones_single_mem", old$build, envir = globalenv())
    assign(".assign_teams_core_mem", old$assign, envir = globalenv())
    pkg$cache_dir <- old$dirs$cache_dir
    pkg$assign_cache_dir <- old$dirs$assign_cache_dir
    unlink(tmp, recursive = TRUE)
  }, envir = testthat::teardown_env())
})

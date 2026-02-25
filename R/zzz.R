#' @keywords internal
"_PACKAGE"

# Suppress R CMD check NOTEs for dplyr/tidyr/ompr/ggplot2 column references
utils::globalVariables(c(
  "origin_id", "destination_id", "distance",
  "i", "j",
  "x", "y", "p",
  "center_id", "tract_id", "distance_to_center",
  "zone_id", "partition_id",
  "xend", "yend",
  "data",  # tidyr::nest() creates this variable
  "surveyzones_build_zones_mem"  # memoised version created in .onLoad
))

# Package-level environment for mutable state
.pkg <- new.env(parent = emptyenv())

.onLoad <- function(libname, pkgname) {
  # ROI requires plugins to be explicitly loaded for solver registration
  requireNamespace("ROI.plugin.glpk", quietly = TRUE)

  # Register solver-specific control mappings for ROI (cf. orce package)
  try(ROI::ROI_plugin_register_solver_control(
    "cbc", "ratio", "rel_tol"
  ), silent = TRUE)
  try(ROI::ROI_plugin_register_solver_control(
    "highs", "mip_rel_gap", "rel_tol"
  ), silent = TRUE)

  # Disk-backed memoisation for the zone solver
  cache_dir <- file.path(tools::R_user_dir("surveyzones", which = "cache"), "build_zones")
  dir.create(cache_dir, showWarnings = FALSE, recursive = TRUE)
  .pkg$cache_dir <- cache_dir

  surveyzones_build_zones_mem <<- memoise::memoise(
    .surveyzones_build_zones_impl,
    cache = cachem::cache_disk(dir = cache_dir)
  )
}

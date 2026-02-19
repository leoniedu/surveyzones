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
  "data"  # tidyr::nest() creates this variable
))

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
}

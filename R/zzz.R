#' @import data.table
#' @keywords internal
"_PACKAGE"

# Ensure data.table NSE works inside the package
.datatable.aware <- TRUE

# Suppress R CMD check NOTEs for data.table/ompr/ggplot2 column references
utils::globalVariables(c(
  "origin_id", "destination_id", "travel_time",
  "i", "j",
  "x", "y", "p",
  "center_id", "tract_id", "travel_time_to_center",
  "zone_id", "partition_id",
  "xend", "yend"
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

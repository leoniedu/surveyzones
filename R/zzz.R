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
}

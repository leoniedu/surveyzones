#' Create a surveyzones_zone_groups Object
#'
#' Constructor for the S3 class returned by [surveyzones_group_zones()].
#'
#' @param zones A tibble with columns `zone_id`, `partition_id`,
#'   `center_tract_id`, `n_tracts`, `group_id`, `ggg`, `zzz`, `label`.
#' @param groups A tibble with columns `group_id`, `partition_id`,
#'   `center_zone_id`, `n_zones`, `diameter`, `ggg`.
#' @param sequence A tibble with columns `zone_id`, `tract_id`,
#'   `visit_order`.
#' @param parameters A named list of parameters (`D_max`, `solver`,
#'   `timestamp`).
#' @param diagnostics A named list of solver diagnostics.
#'
#' @return A `surveyzones_zone_groups` object.
#' @keywords internal
new_surveyzones_zone_groups <- function(zones, groups, sequence,
                                        parameters, diagnostics) {
  structure(
    list(
      zones = zones,
      groups = groups,
      sequence = sequence,
      parameters = parameters,
      diagnostics = diagnostics
    ),
    class = "surveyzones_zone_groups"
  )
}


#' @export
print.surveyzones_zone_groups <- function(x, ...) {
  n_partitions <- length(unique(x$zones$partition_id))
  n_zones <- nrow(x$zones)
  n_groups <- nrow(x$groups)

  cli::cli_h1("surveyzones zone groups")
  cli::cli_bullets(c(
    "*" = "{n_partitions} partition{?s}, {n_zones} zone{?s}, {n_groups} group{?s}",
    "*" = "D_max = {x$parameters$D_max}"
  ))

  if (n_groups > 0) {
    dm <- x$groups$diameter
    cli::cli_bullets(c(
      "*" = "Group diameter: {round(min(dm, na.rm = TRUE), 2)}–{round(max(dm, na.rm = TRUE), 2)} (mean {round(mean(dm, na.rm = TRUE), 2)})"
    ))
  }

  invisible(x)
}

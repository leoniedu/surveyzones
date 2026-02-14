#' Plot Zone Assignments
#'
#' Visualise zone assignments on a map.  Requires tract geometries
#' to be provided as an sf object.
#'
#' @param plan A `surveyzones_plan` object.
#' @param tracts_sf An sf object with a `tract_id` column and polygon
#'   (or point) geometries.
#' @param show_centers Logical.  Highlight zone centers? Default `TRUE`.
#' @param show_sequence Logical.  Draw visit-order arrows?
#'   Requires that [surveyzones_sequence()] has been called on `plan`.
#'   Default `FALSE`.
#'
#' @return A `ggplot` object.
#'
#' @export
surveyzones_plot_zones <- function(plan, tracts_sf,
                                   show_centers = TRUE,
                                   show_sequence = FALSE) {
  if (!inherits(plan, "surveyzones_plan")) {
    cli::cli_abort("{.arg plan} must be a {.cls surveyzones_plan} object.")
  }
  if (!inherits(tracts_sf, "sf")) {
    cli::cli_abort("{.arg tracts_sf} must be an sf object.")
  }
  if (!"tract_id" %in% names(tracts_sf)) {
    cli::cli_abort("{.arg tracts_sf} must contain a {.field tract_id} column.")
  }

  # Join assignments to geometries
  plot_data <- merge(
    tracts_sf,
    plan$assignments[, c("tract_id", "zone_id", "partition_id")],
    by = "tract_id"
  )

  p <- ggplot2::ggplot(plot_data) +
    ggplot2::geom_sf(ggplot2::aes(fill = zone_id), alpha = 0.6) +
    ggplot2::theme_minimal() +
    ggplot2::labs(fill = "Zone") +
    ggplot2::theme(legend.position = "bottom")

  # Centers
  if (show_centers) {
    centers <- plan$assignments[
      plan$assignments$tract_id == plan$assignments$center_id,
    ]
    center_geom <- merge(tracts_sf, centers[, "tract_id"], by = "tract_id")

    if (nrow(center_geom) > 0) {
      center_points <- sf::st_centroid(center_geom)
      p <- p + ggplot2::geom_sf(
        data = center_points,
        shape = 18, size = 3, colour = "black"
      )
    }
  }

  # Sequence arrows
  if (show_sequence) {
    if (is.null(plan$sequence)) {
      cli::cli_warn("No sequence found; run {.fun surveyzones_sequence} first.")
    } else {
      seq_data <- merge(
        plan$sequence,
        sf::st_centroid(tracts_sf)[, c("tract_id", "geometry")],
        by = "tract_id"
      )
      seq_data <- seq_data[order(seq_data$zone_id, seq_data$visit_order), ]
      coords <- sf::st_coordinates(sf::st_as_sf(seq_data))
      seq_data$x <- coords[, 1]
      seq_data$y <- coords[, 2]

      # Build segments
      zones_seq <- split(seq_data, seq_data$zone_id)
      segments <- lapply(zones_seq, function(z) {
        if (nrow(z) < 2) return(NULL)
        data.frame(
          x = z$x[-nrow(z)],
          y = z$y[-nrow(z)],
          xend = z$x[-1],
          yend = z$y[-1],
          zone_id = z$zone_id[-nrow(z)]
        )
      })
      segments <- do.call(rbind, segments)

      if (!is.null(segments) && nrow(segments) > 0) {
        p <- p + ggplot2::geom_segment(
          data = segments,
          ggplot2::aes(x = x, y = y, xend = xend, yend = yend),
          arrow = ggplot2::arrow(length = ggplot2::unit(0.15, "cm")),
          linewidth = 0.5, colour = "grey30"
        )
      }
    }
  }

  # Facet by partition
  n_parts <- length(unique(plot_data$partition_id))
  if (n_parts > 1) {
    p <- p + ggplot2::facet_wrap(~partition_id)
  }

  p
}

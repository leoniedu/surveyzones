utils::globalVariables(c("total_workload", "diameter", "n_tracts", "metric", "value"))

#' Plot Connectivity Graph
#'
#' Visualise the distance connectivity graph (edges within D_max).
#' Shows how tracts are connected via the routing network.
#'
#' @param tracts_sf An sf object with POINT or POLYGON geometries and `tract_id`.
#' @param sparse_distances A sparse distances tibble with `origin_id`,
#'   `destination_id`, and `distance` columns.
#' @param D_max Numeric scalar. Only show edges with distance <= D_max.
#'
#' @return A `ggplot` object.
#' @export
surveyzones_plot_connectivity <- function(tracts_sf, sparse_distances, D_max) {
  if (!inherits(tracts_sf, "sf")) {
    cli::cli_abort("{.arg tracts_sf} must be an sf object.")
  }

  # Filter distances by D_max
  edges <- sparse_distances |>
    dplyr::filter(distance <= D_max)

  if (nrow(edges) == 0) {
    cli::cli_abort("No distance pairs found within D_max = {D_max}")
  }

  # Convert to line geometries
  tract_coords <- sf::st_centroid(tracts_sf) |>
    sf::st_coordinates()
  rownames(tract_coords) <- tracts_sf$tract_id

  edge_lines <- lapply(seq_len(nrow(edges)), function(i) {
    o_id <- edges$origin_id[i]
    d_id <- edges$destination_id[i]

    if (o_id %in% rownames(tract_coords) && d_id %in% rownames(tract_coords)) {
      coords <- rbind(
        tract_coords[o_id, ],
        tract_coords[d_id, ]
      )
      sf::st_linestring(as.matrix(coords))
    } else {
      NULL
    }
  })

  edge_lines <- edge_lines[!sapply(edge_lines, is.null)]

  if (length(edge_lines) == 0) {
    cli::cli_abort("Could not create edge geometries.")
  }

  edge_sfc <- sf::st_sfc(edge_lines, crs = sf::st_crs(tracts_sf))
  edge_sf <- sf::st_sf(geometry = edge_sfc)

  # Plot
  p <- ggplot2::ggplot() +
    ggplot2::geom_sf(data = edge_sf, colour = "grey50", linewidth = 0.3, alpha = 0.5) +
    ggplot2::geom_sf(
      data = sf::st_centroid(tracts_sf),
      size = 2, colour = "steelblue"
    ) +
    ggplot2::theme_minimal() +
    ggplot2::labs(
      title = sprintf("Connectivity Graph (D_max = %g)", D_max),
      caption = sprintf("%d tracts, %d edges", nrow(tracts_sf), nrow(edges))
    )

  p
}


#' Plot Zone Statistics
#'
#' Create summary statistics visualizations for zones.
#'
#' @param plan A `surveyzones_plan` object.
#' @param type Character. Type of plot: `"workload"` (default), `"diameter"`,
#'   `"size"`, or `"all"` for a faceted summary.
#'
#' @return A `ggplot` object.
#' @export
surveyzones_plot_statistics <- function(plan, type = "workload") {
  if (!inherits(plan, "surveyzones_plan")) {
    cli::cli_abort("{.arg plan} must be a {.cls surveyzones_plan} object.")
  }

  type <- match.arg(type, c("workload", "diameter", "size", "all"))
  zones <- plan$zones

  if (type == "workload") {
    ggplot2::ggplot(zones, ggplot2::aes(x = reorder(zone_id, total_workload), y = total_workload)) +
      ggplot2::geom_col(fill = "steelblue", alpha = 0.7) +
      ggplot2::coord_flip() +
      ggplot2::labs(
        title = "Workload per Zone",
        x = "Zone", y = "Total Workload"
      ) +
      ggplot2::theme_minimal()
  } else if (type == "diameter") {
    ggplot2::ggplot(zones, ggplot2::aes(x = reorder(zone_id, diameter), y = diameter)) +
      ggplot2::geom_col(fill = "coral", alpha = 0.7) +
      ggplot2::coord_flip() +
      ggplot2::labs(
        title = "Zone Diameter",
        x = "Zone", y = "Max Diameter (distance units)"
      ) +
      ggplot2::theme_minimal()
  } else if (type == "size") {
    ggplot2::ggplot(zones, ggplot2::aes(x = reorder(zone_id, n_tracts), y = n_tracts)) +
      ggplot2::geom_col(fill = "seagreen", alpha = 0.7) +
      ggplot2::coord_flip() +
      ggplot2::labs(
        title = "Number of Tracts per Zone",
        x = "Zone", y = "Number of Tracts"
      ) +
      ggplot2::theme_minimal()
  } else {
    # Faceted summary
    zones_long <- zones |>
      dplyr::select(zone_id, workload = total_workload, diameter, size = n_tracts) |>
      tidyr::pivot_longer(-zone_id, names_to = "metric", values_to = "value")

    ggplot2::ggplot(zones_long, ggplot2::aes(x = reorder(zone_id, value), y = value, fill = metric)) +
      ggplot2::geom_col(alpha = 0.7) +
      ggplot2::facet_wrap(~metric, scales = "free") +
      ggplot2::coord_flip() +
      ggplot2::labs(title = "Zone Statistics Summary", x = "Zone", y = "Value") +
      ggplot2::theme_minimal() +
      ggplot2::theme(legend.position = "none")
  }
}

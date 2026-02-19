utils::globalVariables(c("total_workload", "diameter", "n_tracts", "metric", "value"))

#' Plot Zone Assignments
#'
#' Visualise zone assignments on a map with optional interactive mode.
#' Requires tract geometries to be provided as an sf object.
#'
#' @importFrom stats reorder
#'
#' @param plan A `surveyzones_plan` object.
#' @param tracts_sf An sf object with a `tract_id` column and polygon
#'   (or point) geometries.
#' @param show_centers Logical.  Highlight zone centers? Default `TRUE`.
#' @param show_sequence Logical.  Draw visit-order arrows?
#'   Requires that [surveyzones_sequence()] has been called on `plan`.
#'   Default `FALSE`.
#' @param interactive Logical.  Deprecated. Use `plot_engine` instead.
#'   Default `FALSE`.
#' @param plot_engine Character.  Engine for rendering: `"ggplot2"` (static),
#'   `"plotly"` (interactive, ggplot-based), `"leaflet"` (interactive, tile-based),
#'   `"mapgl"` (interactive, modern vector tiles), or `"auto"` (default).
#'   When `plot_engine="auto"`, prioritizes mapgl if available, then plotly,
#'   then leaflet, then falls back to ggplot2.
#' @param show_connectivity Logical.  Draw edges showing connectivity graph
#'   within D_max? Requires `sparse_distances` from plan. Default `FALSE`.
#' @param show_labels Logical.  Show zone IDs as text labels? Default `FALSE`.
#'
#' @return A `ggplot` object (ggplot2 engine), `plotly` object (plotly engine),
#'   `leaflet` map object (leaflet engine), or `mapgl` map object (mapgl engine).
#'
#' @export
surveyzones_plot_zones <- function(plan, tracts_sf,
                                   show_centers = TRUE,
                                   show_sequence = FALSE,
                                   interactive = FALSE,
                                   plot_engine = "auto",
                                   show_connectivity = FALSE,
                                   show_labels = FALSE) {
  if (!inherits(plan, "surveyzones_plan")) {
    cli::cli_abort("{.arg plan} must be a {.cls surveyzones_plan} object.")
  }
  if (!inherits(tracts_sf, "sf")) {
    cli::cli_abort("{.arg tracts_sf} must be an sf object.")
  }
  if (!"tract_id" %in% names(tracts_sf)) {
    cli::cli_abort("{.arg tracts_sf} must contain a {.field tract_id} column.")
  }

  # Handle deprecated 'interactive' argument
  if (interactive && plot_engine == "auto") {
    plot_engine <- "plotly"
  }

  # Validate plot_engine argument
  valid_engines <- c("ggplot2", "plotly", "leaflet", "mapgl", "auto")
  if (!plot_engine %in% valid_engines) {
    cli::cli_abort(
      "{.arg plot_engine} must be one of {.val {valid_engines}}, not {.val {plot_engine}}."
    )
  }

  # Resolve "auto" engine choice (prioritize mapgl as modern default)
  if (plot_engine == "auto") {
    if (requireNamespace("mapgl", quietly = TRUE)) {
      plot_engine <- "mapgl"
    } else if (requireNamespace("plotly", quietly = TRUE)) {
      plot_engine <- "plotly"
    } else if (requireNamespace("leaflet", quietly = TRUE)) {
      plot_engine <- "leaflet"
    } else {
      plot_engine <- "ggplot2"
    }
  }

  # Validate requested engine is available
  if (plot_engine == "mapgl" && !requireNamespace("mapgl", quietly = TRUE)) {
    cli::cli_abort(
      "Plot engine {.val mapgl} requires package {.pkg mapgl}."
    )
  }
  if (plot_engine == "plotly" && !requireNamespace("plotly", quietly = TRUE)) {
    cli::cli_abort(
      "Plot engine {.val plotly} requires package {.pkg plotly}. Install with {.code install.packages('plotly')}."
    )
  }
  if (plot_engine == "leaflet" && !requireNamespace("leaflet", quietly = TRUE)) {
    cli::cli_abort(
      "Plot engine {.val leaflet} requires package {.pkg leaflet}. Install with {.code install.packages('leaflet')}."
    )
  }

  # Join assignments to geometries
  plot_data <- merge(
    tracts_sf,
    plan$assignments[, c("tract_id", "zone_id", "partition_id", "center_id")],
    by = "tract_id"
  )

  # For mapgl engine, use dedicated implementation
  if (plot_engine == "mapgl") {
    return(.plot_zones_mapgl(
      plan, plot_data, show_centers, show_sequence, show_labels
    ))
  }

  # For leaflet engine, use dedicated implementation
  if (plot_engine == "leaflet") {
    return(.plot_zones_leaflet(
      plan, plot_data, show_centers, show_sequence, show_labels
    ))
  }

  # Build ggplot visualization (used as base for plotly, or standalone for ggplot2)
  p <- ggplot2::ggplot(plot_data) +
    ggplot2::geom_sf(ggplot2::aes(fill = zone_id), alpha = 0.65, linewidth = 0.2) +
    ggplot2::theme_minimal() +
    ggplot2::labs(fill = "Zone", title = "Zone Assignments") +
    ggplot2::theme(
      legend.position = "bottom",
      plot.title = ggplot2::element_text(face = "bold", size = 12)
    )

  # Zone centers (visual hierarchy: black diamonds)
  if (show_centers) {
    centers <- plan$zones[, c("zone_id", "center_tract_id")]
    center_geom <- merge(
      sf::st_centroid(tracts_sf),
      centers,
      by.x = "tract_id", by.y = "center_tract_id",
      all.x = FALSE
    )

    if (nrow(center_geom) > 0) {
      p <- p + ggplot2::geom_sf(
        data = center_geom,
        shape = 18, size = 4, colour = "black", stroke = 0.5
      )
    }
  }

  # Zone labels (optional)
  if (show_labels) {
    zone_labels <- plan$zones[, c("zone_id")]
    zone_geom <- merge(
      sf::st_centroid(tracts_sf),
      plan$assignments[, c("tract_id", "zone_id")],
      by = "tract_id"
    )
    zone_geom <- unique(zone_geom[, c("zone_id", "geometry")])

    if (nrow(zone_geom) > 0) {
      p <- p + ggplot2::geom_sf_text(
        data = zone_geom,
        ggplot2::aes(label = zone_id),
        size = 3, colour = "white", fontface = "bold"
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

  # Facet by partition (only for ggplot2 static plots)
  if (plot_engine == "ggplot2") {
    n_parts <- length(unique(plot_data$partition_id))
    if (n_parts > 1) {
      p <- p + ggplot2::facet_wrap(~partition_id)
    }
    return(p)
  }

  # Convert to plotly if requested (must be done after plot is fully constructed)
  if (plot_engine == "plotly") {
    return(.plot_zones_plotly(p, plot_data))
  }

  p
}


#' Interactive Zone Visualization with Plotly
#'
#' Convert static ggplot to interactive plotly visualization with hover details.
#' @keywords internal
.plot_zones_plotly <- function(static_plot, plot_data) {
  plotly::ggplotly(static_plot, tooltip = c("fill", "x", "y")) |>
    plotly::layout(
      hovermode = "closest",
      title = list(text = "Interactive Zone Map (hover for details)")
    )
}

#' Interactive Zone Visualization with Leaflet (Fallback)
#'
#' Render zones on an interactive leaflet map (used if maplibre not available).
#' @keywords internal
.plot_zones_leaflet <- function(plan, plot_data, show_centers,
                                show_sequence, show_labels) {
  m <- leaflet::leaflet(plot_data) |>
    leaflet::addTiles(urlTemplate = leaflet::providers$OpenStreetMap.Mapnik) |>
    leaflet::addPolygons(
      fillColor = ~leaflet::colorFactor(
        palette = "Set1",
        domain = plot_data$zone_id
      )(zone_id),
      fillOpacity = 0.6,
      weight = 1,
      color = "black",
      popup = ~paste0(
        "<b>Zone:</b> ", zone_id, "<br/>",
        "<b>Tract:</b> ", tract_id
      )
    )

  # Centers as markers
  if (show_centers) {
    center_coords <- plot_data |>
      subset(!duplicated(zone_id)) |>
      subset(tract_id %in% plan$zones$center_tract_id)

    if (nrow(center_coords) > 0) {
      center_pts <- sf::st_centroid(center_coords)
      coords <- sf::st_coordinates(center_pts)

      m <- m |>
        leaflet::addMarkers(
          lng = coords[, 1],
          lat = coords[, 2],
          popup = ~paste0("<b>Zone Center:</b> ", zone_id),
          icon = leaflet::makeIcon(
            iconUrl = "https://raw.githubusercontent.com/pointhi/leaflet-color-markers/master/img/marker-icon-2x-black.png",
            iconWidth = 25, iconHeight = 41
          )
        )
    }
  }

  # Add summary control
  m <- m |>
    leaflet::addControl(
      html = htmltools::HTML(
        sprintf(
          "<b>Zone Summary</b><br/>Zones: %d<br/>Tracts: %d",
          length(unique(plot_data$zone_id)),
          nrow(plot_data)
        )
      ),
      position = "topright"
    )

  m
}


#' Interactive Zone Visualization with MapGL
#'
#' Render zones on an interactive maplibre map using modern vector tiles.
#' @keywords internal
.plot_zones_mapgl <- function(plan, plot_data, show_centers,
                              show_sequence, show_labels) {
  # Get bounds from data for map initialization
  bounds <- sf::st_bbox(plot_data)

  # Initialize maplibre map
  m <- mapgl::maplibre(bounds = plot_data) |>
    mapgl::add_fill_layer(
      id = "zones",
      source = plot_data,
      fill_color = "zone_id",
      fill_opacity = 0.6,
      stroke_color = "#000000",
      stroke_width = 1
    )

  # Add zone labels if requested
  if (show_labels) {
    # Get zone centroids with zone_id
    zone_labels <- plot_data |>
      dplyr::group_by(zone_id) |>
      dplyr::slice(1) |>
      dplyr::ungroup()

    zone_centroids <- sf::st_centroid(zone_labels)

    m <- m |>
      mapgl::add_symbol_layer(
        id = "zone-labels",
        source = zone_centroids,
        text_field = "zone_id",
        text_size = 12,
        text_color = "#ffffff"
      )
  }

  # Add zone centers as circle markers if requested
  if (show_centers) {
    center_coords <- plot_data |>
      dplyr::group_by(zone_id) |>
      dplyr::slice(1) |>
      dplyr::ungroup() |>
      dplyr::filter(tract_id %in% plan$zones$center_tract_id)

    if (nrow(center_coords) > 0) {
      center_pts <- sf::st_centroid(center_coords)

      m <- m |>
        mapgl::add_circle_layer(
          id = "zone-centers",
          source = center_pts,
          circle_color = "#000000",
          circle_radius = 6,
          circle_stroke_width = 2,
          circle_stroke_color = "#ffffff"
        )
    }
  }

  m
}


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

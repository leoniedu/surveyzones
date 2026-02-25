#' Compute representative points for census tracts
#'
#' Given household address data with coordinates, computes one representative
#' point per census tract using kernel density estimation. This is the main
#' entry point for generating access points compatible with
#' [surveyzones_compute_sparse_distances()].
#'
#' @param households A data frame (or tibble) with household address data. Must
#'   contain columns for tract identifier, latitude, and longitude. Typically
#'   the output of `cnefetools::read_cnefe()` (collected) or any data frame
#'   with geocoded addresses.
#' @param tract_col Character. Column name for the census tract identifier.
#'   Default `"COD_SETOR"` (cnefetools convention).
#' @param lat_col Character. Column name for latitude. Default `"LATITUDE"`.
#' @param lon_col Character. Column name for longitude. Default `"LONGITUDE"`.
#' @param crs Integer. Coordinate reference system for the output.
#'   Default `4674L` (SIRGAS 2000).
#'
#' @return An sf POINT object (CRS as specified) with one row per tract that
#'   has at least one geocoded address. Columns:
#'   \describe{
#'     \item{tract_id}{Census tract identifier (character)}
#'     \item{n_addresses}{Total household addresses in the tract}
#'     \item{n_geocoded}{Number of addresses with valid coordinates}
#'     \item{geometry}{sf POINT geometry — highest-density household location}
#'   }
#'
#' @details
#' For each tract, the function finds the location with the highest
#' concentration of geocoded addresses, using kernel density estimation via
#' [surveyzones_density_point()]. This is more robust than a simple centroid
#' for tracts with irregular shapes or address clusters.
#'
#' Requires the `spatstat.geom` and `spatstat.explore` packages. Install with
#' `install.packages(c("spatstat.geom", "spatstat.explore"))`.
#'
#' @seealso [surveyzones_density_point()], [surveyzones_centroid_point()]
#'
#' @examples
#' \dontrun{
#' library(cnefetools)
#' library(dplyr)
#'
#' # Download CNEFE for Rio de Janeiro city
#' cnefe_rj <- read_cnefe(code_muni = 3304557) |> collect()
#'
#' # Compute one representative point per tract
#' access_pts <- surveyzones_representative_points(cnefe_rj)
#' }
#'
#' @export
surveyzones_representative_points <- function(
    households,
    tract_col = "COD_SETOR",
    lat_col   = "LATITUDE",
    lon_col   = "LONGITUDE",
    crs       = 4674L
) {
  # Validation
  if (!is.data.frame(households)) {
    cli::cli_abort("{.arg households} must be a data frame.")
  }
  missing_cols <- setdiff(c(tract_col, lat_col, lon_col), names(households))
  if (length(missing_cols) > 0) {
    cli::cli_abort(
      "Column{?s} not found in {.arg households}: {.field {missing_cols}}"
    )
  }
  if (nrow(households) == 0) {
    cli::cli_abort("{.arg households} is empty.")
  }

  # Coerce coordinates to numeric
  clean <- households |>
    dplyr::mutate(
      .lat   = suppressWarnings(as.numeric(.data[[lat_col]])),
      .lon   = suppressWarnings(as.numeric(.data[[lon_col]])),
      .tract = as.character(.data[[tract_col]])
    )

  # Summarise per tract
  tract_summary <- clean |>
    dplyr::summarise(
      n_addresses = dplyr::n(),
      n_geocoded  = sum(!is.na(.data$.lat) & !is.na(.data$.lon)),
      .by = ".tract"
    )

  tracts_with_coords <- tract_summary |>
    dplyr::filter(.data$n_geocoded > 0)

  tracts_no_coords <- tract_summary |>
    dplyr::filter(.data$n_geocoded == 0)

  n_total   <- nrow(tract_summary)
  n_with    <- nrow(tracts_with_coords)
  n_without <- nrow(tracts_no_coords)

  cli::cli_h2("Representative Points")
  cli::cli_alert_info("Tracts total: {n_total}")
  cli::cli_alert_info("With coordinates: {n_with}")

  if (n_without > 0) {
    cli::cli_alert_warning(
      "{n_without} tract{?s} have zero geocoded addresses"
    )
  }

  if (n_with == 0) {
    cli::cli_abort("No tracts have geocoded addresses.")
  }

  # Filter to geocoded addresses and convert to sf
  points_sf <- clean |>
    dplyr::filter(
      !is.na(.data$.lat),
      !is.na(.data$.lon),
      .data$.tract %in% tracts_with_coords$.tract
    ) |>
    dplyr::select(".tract", ".lat", ".lon") |>
    dplyr::mutate(n = 1L) |>
    sf::st_as_sf(coords = c(".lon", ".lat"), crs = crs)

  cli::cli_alert_info("Computing density points for {n_with} tract{?s}...")

  repr_points <- surveyzones_density_point(points_sf, .tract)

  result <- repr_points |>
    dplyr::rename(tract_id = ".tract") |>
    dplyr::left_join(
      tracts_with_coords |>
        dplyr::select(
          tract_id = ".tract",
          "n_addresses",
          "n_geocoded"
        ),
      by = dplyr::join_by(tract_id)
    )

  cli::cli_alert_success(
    "Done: {nrow(result)} representative point{?s} computed"
  )

  result
}


#' Find the highest-density point within each spatial unit
#'
#' Uses kernel density estimation (spatstat) to identify the point with the
#' highest concentration of households within each geographic unit (e.g., census
#' tract). This is useful for finding the location where dwellings are most
#' concentrated.
#'
#' @param points An sf POINT object with a column `n` (weight/count) and a
#'   column identified by `geoid`.
#' @param geoid <[`data-masked`][dplyr::dplyr_data_masking]> Column identifying
#'   the geographic unit (e.g., tract code).
#'
#' @return An sf POINT object with columns `geoid` and geometry, one row per
#'   spatial unit, in the same CRS as the input.
#'
#' @details
#' For units with a single point, that point is returned directly without
#' invoking spatstat. For units with multiple points, kernel density estimation
#' is used with a bandwidth of 10\% of the bounding box extent (minimum 30m).
#'
#' Coordinates are projected to SIRGAS 2000 / Brazil Polyconic (EPSG:5880) for
#' the density calculation, then the selected point is returned in the original
#' CRS.
#'
#' Requires `spatstat.geom` and `spatstat.explore` (suggested, not imported).
#' Install with `install.packages(c("spatstat.geom", "spatstat.explore"))`.
#'
#' @seealso [surveyzones_representative_points()], [surveyzones_centroid_point()]
#'
#' @examples
#' \dontrun{
#' library(sf)
#' pts <- st_as_sf(
#'   data.frame(
#'     tract = rep(c("A", "B"), each = 10),
#'     n     = rpois(20, 3),
#'     lon   = c(rnorm(10, -43.2, 0.01), rnorm(10, -43.3, 0.01)),
#'     lat   = c(rnorm(10, -22.9, 0.01), rnorm(10, -23.0, 0.01))
#'   ),
#'   coords = c("lon", "lat"),
#'   crs    = 4674
#' )
#' surveyzones_density_point(pts, tract)
#' }
#'
#' @export
surveyzones_density_point <- function(points, geoid) {
  geoid_nm <- rlang::as_name(rlang::enquo(geoid))

  counts <- points |>
    sf::st_drop_geometry() |>
    dplyr::count(.data[[geoid_nm]], name = ".n_pts")

  single_ids <- counts |>
    dplyr::filter(.data$.n_pts == 1) |>
    dplyr::pull(.data[[geoid_nm]])

  multi_ids <- counts |>
    dplyr::filter(.data$.n_pts > 1) |>
    dplyr::pull(.data[[geoid_nm]])

  result_single <- if (length(single_ids) > 0) {
    points |>
      dplyr::filter(.data[[geoid_nm]] %in% single_ids) |>
      dplyr::select(dplyr::all_of(geoid_nm))
  }

  result_multi <- if (length(multi_ids) > 0) {
    pts_proj <- sf::st_transform(points, crs = 5880L)

    use_progress <- length(multi_ids) > 50
    if (use_progress) {
      bar_id <- cli::cli_progress_bar(
        "Computing density points",
        total = length(multi_ids)
      )
    }

    results <- purrr::map(multi_ids, \(id) {
      subset_proj <- pts_proj |>
        dplyr::filter(.data[[geoid_nm]] == id)

      idx <- .density_point_index(subset_proj)
      if (use_progress) {
        cli::cli_progress_update(id = bar_id)
      }

      points |>
        dplyr::filter(.data[[geoid_nm]] == id) |>
        dplyr::slice(idx) |>
        dplyr::select(dplyr::all_of(geoid_nm))
    })

    if (use_progress) {
      cli::cli_progress_done(id = bar_id)
    }
    dplyr::bind_rows(results)
  }

  dplyr::bind_rows(result_single, result_multi)
}


#' Find the weighted centroid within each spatial unit
#'
#' Computes the weighted centroid of points within each geographic unit using
#' `n` as weights. A simpler, faster alternative to [surveyzones_density_point()]
#' when kernel density estimation is not needed.
#'
#' @inheritParams surveyzones_density_point
#'
#' @return An sf POINT object with columns `geoid` and geometry, one row per
#'   spatial unit, in the same CRS as the input.
#'
#' @seealso [surveyzones_representative_points()], [surveyzones_density_point()]
#'
#' @examples
#' \dontrun{
#' library(sf)
#' pts <- st_as_sf(
#'   data.frame(
#'     tract = rep(c("A", "B"), each = 5),
#'     n     = rpois(10, 3),
#'     lon   = c(rnorm(5, -43.2, 0.01), rnorm(5, -43.3, 0.01)),
#'     lat   = c(rnorm(5, -22.9, 0.01), rnorm(5, -23.0, 0.01))
#'   ),
#'   coords = c("lon", "lat"),
#'   crs    = 4674
#' )
#' surveyzones_centroid_point(pts, tract)
#' }
#'
#' @export
surveyzones_centroid_point <- function(points, geoid) {
  geoid_nm <- rlang::as_name(rlang::enquo(geoid))

  original_crs <- sf::st_crs(points)
  coords <- sf::st_coordinates(points)

  points |>
    sf::st_drop_geometry() |>
    dplyr::mutate(.lon = coords[, 1], .lat = coords[, 2]) |>
    dplyr::summarise(
      .lon = stats::weighted.mean(.data$.lon, w = .data$n),
      .lat = stats::weighted.mean(.data$.lat, w = .data$n),
      .by  = dplyr::all_of(geoid_nm)
    ) |>
    sf::st_as_sf(coords = c(".lon", ".lat"), crs = original_crs)
}


# Internal -----------------------------------------------------------------

#' Find index of highest-density point in a projected sf subset
#'
#' @param subset_proj sf POINT object already projected to a metric CRS
#'   (e.g., EPSG:5880). Must have a column `n` with weights.
#' @return Integer index (length 1) of the highest-density point.
#' @noRd
.density_point_index <- function(subset_proj) {
  rlang::check_installed(
    c("spatstat.geom", "spatstat.explore"),
    reason = "for kernel density estimation in surveyzones_density_point()"
  )

  xy      <- sf::st_coordinates(subset_proj)
  x       <- xy[, 1]
  y       <- xy[, 2]
  weights <- subset_proj$n

  x_range <- range(x)
  y_range <- range(y)

  if (diff(x_range) < 1e-6 && diff(y_range) < 1e-6) {
    return(which.max(weights))
  }

  sigma <- pmax(diff(x_range) / 10, diff(y_range) / 10, 30)

  ppp_obj <- spatstat.geom::ppp(
    x      = x,
    y      = y,
    window = spatstat.geom::owin(xrange = x_range, yrange = y_range),
    check  = FALSE
  )

  dp <- spatstat.explore::density.ppp(
    ppp_obj,
    sigma   = sigma,
    weights = weights,
    at      = "points"
  )

  if (is.list(dp)) {
    dp <- as.numeric(dp)
  }

  which.max(dp)
}

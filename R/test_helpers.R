#' Create Synthetic Tract Data for Testing
#'
#' Generates a small set of census tracts with point geometries
#' arranged in a grid pattern, mimicking a survey sample structure.
#' Column names follow the surveyzones schema (English).
#'
#' @param n_tracts Integer.  Number of tracts to generate.
#' @param n_partitions Integer.  Number of jurisdictions.
#'   Tracts are distributed evenly across partitions.
#' @param bbox Numeric vector of length 4:
#'   `c(xmin, ymin, xmax, ymax)` in longitude/latitude.
#'   Default covers a small area around Salvador, BA.
#' @param seed Integer.  Random seed for reproducibility.
#'
#' @return An sf object with columns `tract_id`,
#'   `expected_service_time`, `partition_id`, and POINT geometry.
#'
#' @examples
#' tracts <- surveyzones_example_tracts(n_tracts = 20)
#' plot(tracts["partition_id"])
#'
#' @export
surveyzones_example_tracts <- function(
    n_tracts = 24L,
    n_partitions = 2L,
    bbox = c(-38.55, -13.05, -38.45, -12.95),
    seed = 42L) {

  set.seed(seed)

  lon <- stats::runif(n_tracts, bbox[1], bbox[3])
  lat <- stats::runif(n_tracts, bbox[2], bbox[4])

  tracts <- tibble::tibble(
    tract_id = sprintf("T%03d", seq_len(n_tracts)),
    expected_service_time = round(stats::runif(n_tracts, 0.5, 2.0), 1),
    partition_id = sprintf("AG%02d", rep_len(seq_len(n_partitions), n_tracts))
  )

  geom <- sf::st_sfc(
    lapply(seq_len(n_tracts), function(i) sf::st_point(c(lon[i], lat[i]))),
    crs = 4326
  )

  sf::st_sf(tracts, geometry = geom)
}

# Tests for surveyzones_density_point() and surveyzones_centroid_point()

make_pts <- function(n_per_group = 10, crs = 4674) {
  set.seed(42)
  sf::st_as_sf(
    data.frame(
      tract = rep(c("A", "B"), each = n_per_group),
      n     = rpois(2 * n_per_group, lambda = 3) + 1L,
      lon   = c(
        rnorm(n_per_group, -38.50, 0.005),
        rnorm(n_per_group, -38.60, 0.005)
      ),
      lat   = c(
        rnorm(n_per_group, -12.90, 0.005),
        rnorm(n_per_group, -13.00, 0.005)
      )
    ),
    coords = c("lon", "lat"),
    crs    = crs
  )
}

# surveyzones_centroid_point -----------------------------------------------

test_that("surveyzones_centroid_point returns sf POINT with correct columns", {
  pts    <- make_pts()
  result <- surveyzones_centroid_point(pts, tract)

  expect_s3_class(result, "sf")
  expect_true("tract" %in% names(result))
  expect_equal(nrow(result), 2)
  expect_equal(sf::st_crs(result), sf::st_crs(pts))
  expect_true(all(sf::st_geometry_type(result) == "POINT"))
})

test_that("surveyzones_centroid_point returns one row per group", {
  pts    <- make_pts(n_per_group = 20)
  result <- surveyzones_centroid_point(pts, tract)
  expect_equal(nrow(result), 2)
  expect_equal(sort(result$tract), c("A", "B"))
})

test_that("surveyzones_centroid_point computes weighted mean coordinates", {
  pts <- sf::st_as_sf(
    data.frame(grp = c("X", "X"), n = c(1, 3), lon = c(0, 10), lat = c(0, 0)),
    coords = c("lon", "lat"), crs = 4326
  )
  result <- surveyzones_centroid_point(pts, grp)
  coords <- sf::st_coordinates(result)
  expect_equal(coords[1, "X"], 7.5, ignore_attr = TRUE)
  expect_equal(coords[1, "Y"], 0,   ignore_attr = TRUE)
})

test_that("surveyzones_centroid_point handles single-point groups", {
  pts <- sf::st_as_sf(
    data.frame(grp = "A", n = 5, lon = -38.5, lat = -12.9),
    coords = c("lon", "lat"), crs = 4674
  )
  result <- surveyzones_centroid_point(pts, grp)
  expect_equal(nrow(result), 1)
  coords <- sf::st_coordinates(result)
  expect_equal(coords[1, "X"], -38.5, ignore_attr = TRUE)
  expect_equal(coords[1, "Y"], -12.9, ignore_attr = TRUE)
})


# surveyzones_density_point (requires spatstat) ----------------------------

test_that("surveyzones_density_point returns sf POINT with correct columns", {
  skip_if_not_installed("spatstat.geom")
  skip_if_not_installed("spatstat.explore")

  pts    <- make_pts()
  result <- surveyzones_density_point(pts, tract)

  expect_s3_class(result, "sf")
  expect_true("tract" %in% names(result))
  expect_equal(nrow(result), 2)
  expect_equal(sf::st_crs(result), sf::st_crs(pts))
  expect_true(all(sf::st_geometry_type(result) == "POINT"))
})

test_that("surveyzones_density_point returns one row per group", {
  skip_if_not_installed("spatstat.geom")
  skip_if_not_installed("spatstat.explore")

  pts    <- make_pts(n_per_group = 20)
  result <- surveyzones_density_point(pts, tract)
  expect_equal(nrow(result), 2)
  expect_equal(sort(result$tract), c("A", "B"))
})

test_that("surveyzones_density_point handles single-point groups without spatstat", {
  skip_if_not_installed("spatstat.geom")
  skip_if_not_installed("spatstat.explore")

  pts <- sf::st_as_sf(
    data.frame(
      grp = c("single", "multi", "multi", "multi"),
      n   = c(5, 1, 10, 1),
      lon = c(-38.5, -38.6, -38.6001, -38.5999),
      lat = c(-12.9, -13.0, -13.0001, -12.9999)
    ),
    coords = c("lon", "lat"), crs = 4674
  )
  result <- surveyzones_density_point(pts, grp)
  expect_equal(nrow(result), 2)
  expect_equal(sort(result$grp), c("multi", "single"))

  single_result <- result[result$grp == "single", ]
  single_input  <- pts[pts$grp == "single", ]
  expect_equal(sf::st_coordinates(single_result), sf::st_coordinates(single_input))
})

test_that("surveyzones_density_point selects high-weight cluster", {
  skip_if_not_installed("spatstat.geom")
  skip_if_not_installed("spatstat.explore")

  set.seed(123)
  pts <- sf::st_as_sf(
    data.frame(
      grp = rep("X", 11),
      n   = c(rep(10, 10), 1),
      lon = c(rnorm(10, -38.5, 0.001), -39.0),
      lat = c(rnorm(10, -12.9, 0.001), -13.5)
    ),
    coords = c("lon", "lat"), crs = 4674
  )
  result <- surveyzones_density_point(pts, grp)
  coords <- sf::st_coordinates(result)
  expect_true(coords[1, "X"] > -38.6)
  expect_true(coords[1, "Y"] > -13.0)
})

test_that("surveyzones_density_point preserves CRS", {
  skip_if_not_installed("spatstat.geom")
  skip_if_not_installed("spatstat.explore")

  pts    <- make_pts(crs = 4326)
  result <- surveyzones_density_point(pts, tract)
  expect_equal(sf::st_crs(result)$epsg, 4326L)
})

test_that("surveyzones_density_point works with different geoid column names", {
  skip_if_not_installed("spatstat.geom")
  skip_if_not_installed("spatstat.explore")

  pts <- sf::st_as_sf(
    data.frame(
      setor_codigo = rep(c("001", "002"), each = 5),
      n            = rep(1, 10),
      lon          = c(rnorm(5, -38.5, 0.01), rnorm(5, -38.6, 0.01)),
      lat          = c(rnorm(5, -12.9, 0.01), rnorm(5, -13.0, 0.01))
    ),
    coords = c("lon", "lat"), crs = 4674
  )
  result <- surveyzones_density_point(pts, setor_codigo)
  expect_true("setor_codigo" %in% names(result))
  expect_equal(nrow(result), 2)
})

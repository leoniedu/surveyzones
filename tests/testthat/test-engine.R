test_that("haversine engine returns correct dimensions", {
  pts <- sf::st_as_sf(
    data.frame(id = 1:3, lon = c(-38.5, -38.51, -38.52),
               lat = c(-13.0, -13.01, -13.02)),
    coords = c("lon", "lat"), crs = 4326
  )

  engine <- surveyzones_engine_haversine(units = "km")
  mat <- engine(pts[1:2, ], pts)

  expect_true(is.matrix(mat))
  expect_equal(nrow(mat), 2)
  expect_equal(ncol(mat), 3)
  expect_true(all(mat >= 0))
})

test_that("haversine engine respects units", {
  pts <- sf::st_as_sf(
    data.frame(lon = c(-38.5, -38.6), lat = c(-13.0, -13.0)),
    coords = c("lon", "lat"), crs = 4326
  )

  km_engine <- surveyzones_engine_haversine(units = "km")
  m_engine <- surveyzones_engine_haversine(units = "m")

  d_km <- km_engine(pts[1, ], pts[2, ])
  d_m <- m_engine(pts[1, ], pts[2, ])

  expect_equal(d_m[1, 1], d_km[1, 1] * 1000, tolerance = 0.01)
})

test_that("osrm engine factory returns a valid engine", {
  skip_if_not_installed("osrm")

  engine <- surveyzones_engine_osrm(measure = "duration")
  expect_true(is.function(engine))
  expect_silent(validate_engine(engine))
})

test_that("osrm engine returns correct dimensions", {
  skip_if_not_installed("osrm")
  skip_if_offline()
  skip_on_cran()

  pts <- sf::st_as_sf(
    data.frame(id = 1:3, lon = c(-43.17, -43.18, -43.19),
               lat = c(-22.91, -22.92, -22.93)),
    coords = c("lon", "lat"), crs = 4326
  )

  engine <- surveyzones_engine_osrm(measure = "duration")
  mat <- engine(pts[1:2, ], pts)

  expect_true(is.matrix(mat))
  expect_equal(nrow(mat), 2)
  expect_equal(ncol(mat), 3)
  expect_true(all(is.finite(mat)))
})

test_that("osrm engine respects measure parameter", {
  skip_if_not_installed("osrm")
  skip_if_offline()
  skip_on_cran()

  pts <- sf::st_as_sf(
    data.frame(id = 1:2, lon = c(-43.17, -43.18),
               lat = c(-22.91, -22.92)),
    coords = c("lon", "lat"), crs = 4326
  )

  dur_engine <- surveyzones_engine_osrm(measure = "duration")
  dist_engine <- surveyzones_engine_osrm(measure = "distance")

  dur_mat <- dur_engine(pts[1, ], pts[2, ])
  dist_mat <- dist_engine(pts[1, ], pts[2, ])

  # Duration in minutes, distance in meters — they should differ
  expect_false(dur_mat[1, 1] == dist_mat[1, 1])
})

test_that("validate_engine rejects non-functions", {
  expect_error(
    validate_engine("not_a_function"),
    "must be a function"
  )
})

test_that("validate_engine rejects wrong arity", {
  expect_error(
    validate_engine(function(x) x),
    "at least 2 arguments"
  )
})

test_that("custom engine is accepted", {
  my_engine <- function(origins, destinations) {
    matrix(1.0, nrow = nrow(origins), ncol = nrow(destinations))
  }

  expect_silent(validate_engine(my_engine))
})

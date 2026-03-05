test_that("sparse distances have correct schema", {
  pts <- sf::st_as_sf(
    data.frame(
      tract_id = c("A", "B", "C"),
      lon = c(-38.5, -38.51, -38.52),
      lat = c(-13.0, -13.01, -13.02)
    ),
    coords = c("lon", "lat"), crs = 4326
  )

  result <- surveyzones_compute_sparse_distances(
    access_points = pts,
    chunk_size = 2
  )

  expect_s3_class(result, "tbl_df")
  expect_true(all(c("origin_id", "destination_id", "distance") %in% names(result)))
})

test_that("sparse distances exclude self-pairs", {
  pts <- sf::st_as_sf(
    data.frame(
      tract_id = c("A", "B"),
      lon = c(-38.5, -38.51),
      lat = c(-13.0, -13.01)
    ),
    coords = c("lon", "lat"), crs = 4326
  )

  result <- surveyzones_compute_sparse_distances(pts)
  self <- result |> dplyr::filter(origin_id == destination_id)
  expect_equal(nrow(self), 0)
})

test_that("sparse distances can be filtered by D_max", {
  pts <- sf::st_as_sf(
    data.frame(
      tract_id = c("A", "B", "C"),
      lon = c(-38.5, -38.51, -40.0),  # C is far away
      lat = c(-13.0, -13.01, -13.0)
    ),
    coords = c("lon", "lat"), crs = 4326
  )

  all_distances <- surveyzones_compute_sparse_distances(pts)
  # Caller now filters by D_max
  result <- all_distances |> dplyr::filter(distance <= 5)

  expect_true(all(result$distance <= 5))
  # C should not appear in connections to A or B when filtered
  far_pairs <- result |>
    dplyr::filter(
      (origin_id == "C" & destination_id %in% c("A", "B")) |
      (destination_id == "C" & origin_id %in% c("A", "B"))
    )
  expect_equal(nrow(far_pairs), 0)
})

test_that("precomputed distances filter correctly", {
  dt <- data.frame(
    origin_id = c("A", "A", "B"),
    destination_id = c("B", "C", "C"),
    distance = c(3, 10, 5)
  )

  result <- surveyzones_precomputed_distances(dt, D_max = 6)
  expect_equal(nrow(result), 2)  # only A-B (3) and B-C (5)
  expect_true(all(result$distance <= 6))
})

test_that("precomputed distances reject missing columns", {
  dt <- data.frame(from = "A", to = "B", dist = 1)
  expect_error(
    surveyzones_precomputed_distances(dt, D_max = 5),
    "missing column"
  )
})

test_that("validate_access_points rejects non-sf", {
  expect_error(
    validate_access_points(data.frame(tract_id = "A")),
    "must be an sf object"
  )
})

test_that("validate_access_points rejects missing tract_id", {
  pts <- sf::st_as_sf(
    data.frame(id = 1, lon = -38.5, lat = -13.0),
    coords = c("lon", "lat"), crs = 4326
  )
  expect_error(
    validate_access_points(pts),
    "tract_id"
  )
})

test_that("validate_access_points rejects duplicate tract_ids", {
  pts <- sf::st_as_sf(
    data.frame(tract_id = c("A", "A"), lon = c(-38.5, -38.6), lat = c(-13, -13)),
    coords = c("lon", "lat"), crs = 4326
  )
  expect_error(
    validate_access_points(pts),
    "must be unique"
  )
})

test_that("complete_distances fills missing pairs with haversine", {
  pts <- sf::st_as_sf(
    data.frame(
      tract_id = c("A", "B", "C"),
      lon = c(-38.50, -38.51, -38.52),
      lat = c(-13.00, -13.01, -13.02)
    ),
    coords = c("lon", "lat"), crs = 4326
  )

  # Only A-B pair exists; A-C and B-C are missing
  sparse <- tibble::tibble(
    origin_id = c("A", "B"),
    destination_id = c("B", "A"),
    distance = c(10, 10)
  )

  result <- surveyzones_complete_distances(sparse, pts, speed_kmh = 0.1)

  # Original pairs preserved
  ab <- result |> dplyr::filter(origin_id == "A", destination_id == "B")
  expect_equal(ab$distance, 10)

  # Missing pairs filled (all 6 directed pairs for 3 points should exist)
  expect_equal(nrow(result), 6)

  # Filled pairs have positive distance
  ac <- result |> dplyr::filter(origin_id == "A", destination_id == "C")
  expect_true(ac$distance > 0)
})

test_that("complete_distances preserves all existing pairs", {
  pts <- sf::st_as_sf(
    data.frame(
      tract_id = c("A", "B"),
      lon = c(-38.50, -38.51),
      lat = c(-13.00, -13.01)
    ),
    coords = c("lon", "lat"), crs = 4326
  )

  # Already complete
  sparse <- tibble::tibble(
    origin_id = c("A", "B"),
    destination_id = c("B", "A"),
    distance = c(99, 99)
  )

  result <- surveyzones_complete_distances(sparse, pts, speed_kmh = 0.1)

  # No new rows added, original values preserved
  expect_equal(nrow(result), 2)
  expect_true(all(result$distance == 99))
})

test_that("complete_distances respects speed_kmh parameter", {
  pts <- sf::st_as_sf(
    data.frame(
      tract_id = c("A", "B"),
      lon = c(-38.50, -38.51),
      lat = c(-13.00, -13.01)
    ),
    coords = c("lon", "lat"), crs = 4326
  )

  sparse <- tibble::tibble(
    origin_id = character(0),
    destination_id = character(0),
    distance = numeric(0)
  )

  slow <- surveyzones_complete_distances(sparse, pts, speed_kmh = 0.1)
  fast <- surveyzones_complete_distances(sparse, pts, speed_kmh = 1.0)

  # Slower speed -> larger distance values
  slow_d <- slow$distance[slow$origin_id == "A" & slow$destination_id == "B"]
  fast_d <- fast$distance[fast$origin_id == "A" & fast$destination_id == "B"]
  expect_true(slow_d > fast_d)
  expect_equal(slow_d / fast_d, 10, tolerance = 0.01)
})

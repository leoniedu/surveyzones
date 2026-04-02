test_that("spopt assigns all tracts exactly once", {
  skip_if_not_installed("spopt")

  pts <- sf::st_as_sf(
    data.frame(
      tract_id = sprintf("T%d", 1:6),
      lon = c(-38.50, -38.501, -38.502, -38.51, -38.511, -38.512),
      lat = c(-13.00, -13.001, -13.002, -13.01, -13.011, -13.012)
    ),
    coords = c("lon", "lat"), crs = 4326
  )

  tracts <- data.frame(
    tract_id = sprintf("T%d", 1:6),
    expected_service_time = rep(1, 6),
    partition_id = c("A", "A", "A", "B", "B", "B")
  )

  dists <- surveyzones_compute_sparse_distances(pts)

  plan <- surveyzones_build_zones(
    sparse_distances = dists,
    tracts = tracts,
    D_max = 10,
    max_workload_per_zone = 2,
    enforce_partition = TRUE,
    solver = "spopt",
    use_cache = FALSE
  )

  expect_s3_class(plan, "surveyzones_plan")
  expect_equal(sort(plan$assignments$tract_id), sort(tracts$tract_id))
  expect_equal(anyDuplicated(plan$assignments$tract_id), 0)
})

test_that("spopt no zone exceeds max_workload", {
  skip_if_not_installed("spopt")

  pts <- sf::st_as_sf(
    data.frame(
      tract_id = sprintf("T%d", 1:4),
      lon = c(-38.50, -38.501, -38.502, -38.503),
      lat = c(-13.00, -13.001, -13.002, -13.003)
    ),
    coords = c("lon", "lat"), crs = 4326
  )

  tracts <- data.frame(
    tract_id = sprintf("T%d", 1:4),
    expected_service_time = c(1.5, 1.5, 1.5, 1.5)
  )

  dists <- surveyzones_compute_sparse_distances(pts)

  plan <- surveyzones_build_zones(
    sparse_distances = dists,
    tracts = tracts,
    D_max = 10,
    max_workload_per_zone = 3,
    enforce_partition = FALSE,
    solver = "spopt",
    use_cache = FALSE
  )

  expect_true(all(plan$zones$total_workload <= 3))
})

test_that("spopt partition boundaries are respected", {
  skip_if_not_installed("spopt")

  pts <- sf::st_as_sf(
    data.frame(
      tract_id = sprintf("T%d", 1:4),
      lon = c(-38.50, -38.501, -38.51, -38.511),
      lat = c(-13.00, -13.001, -13.01, -13.011)
    ),
    coords = c("lon", "lat"), crs = 4326
  )

  tracts <- data.frame(
    tract_id = sprintf("T%d", 1:4),
    expected_service_time = c(1, 1, 1, 1),
    partition_id = c("A", "A", "B", "B")
  )

  dists <- surveyzones_compute_sparse_distances(pts)

  plan <- surveyzones_build_zones(
    sparse_distances = dists,
    tracts = tracts,
    D_max = 10,
    enforce_partition = TRUE,
    solver = "spopt",
    use_cache = FALSE
  )

  # Each zone should only have tracts from one partition
  zone_parts <- plan$assignments |>
    dplyr::left_join(
      data.frame(tract_id = tracts$tract_id, partition_id = tracts$partition_id),
      by = "tract_id"
    ) |>
    dplyr::summarise(
      n_parts = dplyr::n_distinct(.data$partition_id.y),
      .by = "zone_id"
    )

  expect_true(all(zone_parts$n_parts == 1))
})

test_that("spopt uncapacitated produces valid plan", {
  skip_if_not_installed("spopt")

  pts <- sf::st_as_sf(
    data.frame(
      tract_id = sprintf("T%d", 1:4),
      lon = c(-38.50, -38.501, -38.502, -38.503),
      lat = c(-13.00, -13.001, -13.002, -13.003)
    ),
    coords = c("lon", "lat"), crs = 4326
  )

  tracts <- data.frame(
    tract_id = sprintf("T%d", 1:4),
    expected_service_time = c(1, 1, 1, 1)
  )

  dists <- surveyzones_compute_sparse_distances(pts)

  plan <- surveyzones_build_zones(
    sparse_distances = dists,
    tracts = tracts,
    D_max = 10,
    enforce_partition = FALSE,
    solver = "spopt",
    use_cache = FALSE
  )

  expect_s3_class(plan, "surveyzones_plan")
  # All nearby tracts should be in 1 zone
  expect_equal(nrow(plan$zones), 1L)
})

test_that("spopt plan works with surveyzones_sequence", {
  skip_if_not_installed("spopt")

  pts <- sf::st_as_sf(
    data.frame(
      tract_id = sprintf("T%d", 1:6),
      lon = c(-38.50, -38.501, -38.502, -38.503, -38.504, -38.505),
      lat = c(-13.00, -13.001, -13.002, -13.003, -13.004, -13.005)
    ),
    coords = c("lon", "lat"), crs = 4326
  )

  tracts <- data.frame(
    tract_id = sprintf("T%d", 1:6),
    expected_service_time = rep(1, 6)
  )

  dists <- surveyzones_compute_sparse_distances(pts)

  plan <- surveyzones_build_zones(
    sparse_distances = dists,
    tracts = tracts,
    D_max = 10,
    max_workload_per_zone = 3,
    enforce_partition = FALSE,
    solver = "spopt",
    use_cache = FALSE
  )

  sequenced <- surveyzones_sequence(plan, dists, complete_distances = FALSE)

  expect_s3_class(sequenced, "surveyzones_plan")
  expect_true(!is.null(sequenced$sequence))
  expect_true(!is.null(sequenced$zone_sequence))
})

test_that("spopt plan records solver parameter", {
  skip_if_not_installed("spopt")

  pts <- sf::st_as_sf(
    data.frame(
      tract_id = sprintf("T%d", 1:3),
      lon = c(-38.50, -38.501, -38.502),
      lat = c(-13.00, -13.001, -13.002)
    ),
    coords = c("lon", "lat"), crs = 4326
  )

  tracts <- data.frame(
    tract_id = sprintf("T%d", 1:3),
    expected_service_time = rep(1, 3)
  )

  dists <- surveyzones_compute_sparse_distances(pts)

  plan <- surveyzones_build_zones(
    sparse_distances = dists,
    tracts = tracts,
    D_max = 10,
    enforce_partition = FALSE,
    solver = "spopt",
    use_cache = FALSE
  )

  expect_equal(plan$parameters$solver, "spopt")
})

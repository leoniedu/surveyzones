test_that("build_zones assigns all tracts exactly once", {
  # Small problem: 6 tracts, close together, 2 partitions
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

  dists <- surveyzones_compute_sparse_distances(pts, D_max = 10)

  plan <- surveyzones_build_zones(
    sparse_distances = dists,
    tracts = tracts,
    D_max = 10,
    max_workload_per_zone = 2,
    enforce_partition = TRUE
  )

  expect_s3_class(plan, "surveyzones_plan")
  # All 6 tracts assigned
  expect_equal(sort(plan$assignments$tract_id), sort(tracts$tract_id))
  # No duplicates
  expect_equal(anyDuplicated(plan$assignments$tract_id), 0)
})

test_that("no zone exceeds max_workload", {
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

  dists <- surveyzones_compute_sparse_distances(pts, D_max = 10)

  plan <- surveyzones_build_zones(
    sparse_distances = dists,
    tracts = tracts,
    D_max = 10,
    max_workload_per_zone = 3,
    enforce_partition = FALSE
  )

  expect_true(all(plan$zones$total_workload <= 3))
})

test_that("partition boundaries are respected", {
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

  dists <- surveyzones_compute_sparse_distances(pts, D_max = 10)

  plan <- surveyzones_build_zones(
    sparse_distances = dists,
    tracts = tracts,
    D_max = 10,
    max_workload_per_zone = 2,
    enforce_partition = TRUE
  )

  # Tracts in partition A should only be in zones from partition A
  a_tracts <- plan$assignments[plan$assignments$partition_id == "A", ]
  b_tracts <- plan$assignments[plan$assignments$partition_id == "B", ]

  a_centers <- unique(a_tracts$center_id)
  b_centers <- unique(b_tracts$center_id)

  # Centers from A should not appear in B and vice versa

  expect_equal(length(intersect(a_centers, b_centers)), 0)
})

test_that("validate_tracts catches problems", {
  expect_error(
    validate_tracts(data.frame(x = 1)),
    "missing column"
  )

  expect_error(
    validate_tracts(data.frame(
      tract_id = c("A", "A"),
      expected_service_time = c(1, 1)
    )),
    "must be unique"
  )

  expect_error(
    validate_tracts(data.frame(
      tract_id = c("A", "B"),
      expected_service_time = c(-1, 1)
    )),
    "must be > 0"
  )
})

test_that("print and summary don't error", {
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

  dists <- surveyzones_compute_sparse_distances(pts, D_max = 10)

  plan <- surveyzones_build_zones(
    sparse_distances = dists,
    tracts = tracts,
    D_max = 10,
    max_workload_per_zone = 2,
    enforce_partition = FALSE
  )

  expect_no_error(print(plan))
  expect_no_error(summary(plan))
})

test_that("validate_solver rejects invalid solvers", {
  expect_error(validate_solver("symphony"), "must be one of")
  expect_error(validate_solver(42), "must be a single character")
  expect_error(validate_solver(c("glpk", "highs")), "must be a single character")
})

test_that("validate_solver accepts valid solvers", {
  expect_invisible(validate_solver("glpk"))
})

test_that("solver parameter is recorded in plan", {
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

  dists <- surveyzones_compute_sparse_distances(pts, D_max = 10)

  plan <- surveyzones_build_zones(
    sparse_distances = dists,
    tracts = tracts,
    D_max = 10,
    max_workload_per_zone = 2,
    enforce_partition = FALSE,
    solver = "glpk"
  )

  expect_equal(plan$parameters$solver, "glpk")
})

test_that("highs solver produces valid zones", {
  skip_if_not_installed("ROI.plugin.highs")

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

  dists <- surveyzones_compute_sparse_distances(pts, D_max = 10)

  plan <- surveyzones_build_zones(
    sparse_distances = dists,
    tracts = tracts,
    D_max = 10,
    max_workload_per_zone = 2,
    enforce_partition = FALSE,
    solver = "highs"
  )

  expect_s3_class(plan, "surveyzones_plan")
  expect_equal(sort(plan$assignments$tract_id), sort(tracts$tract_id))
  expect_true(all(plan$zones$total_workload <= 2))
  expect_equal(plan$parameters$solver, "highs")
})

test_that("cbc solver produces valid zones", {
  skip_if_not_installed("ROI.plugin.cbc")

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

  dists <- surveyzones_compute_sparse_distances(pts, D_max = 10)

  plan <- surveyzones_build_zones(
    sparse_distances = dists,
    tracts = tracts,
    D_max = 10,
    max_workload_per_zone = 2,
    enforce_partition = FALSE,
    solver = "cbc"
  )

  expect_s3_class(plan, "surveyzones_plan")
  expect_equal(sort(plan$assignments$tract_id), sort(tracts$tract_id))
  expect_true(all(plan$zones$total_workload <= 2))
  expect_equal(plan$parameters$solver, "cbc")
})

test_that("target_zone_size determines K correctly", {
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

  dists <- surveyzones_compute_sparse_distances(pts, D_max = 10)

  # target_zone_size = 2 -> K = ceil(6/2) = 3 zones
  plan <- surveyzones_build_zones(
    sparse_distances = dists,
    tracts = tracts,
    D_max = 10,
    target_zone_size = 2,
    enforce_partition = FALSE
  )

  expect_s3_class(plan, "surveyzones_plan")
  expect_equal(nrow(plan$zones), 3L)
  expect_equal(sort(plan$assignments$tract_id), sort(tracts$tract_id))
  expect_equal(plan$parameters$target_zone_size, 2)
})

test_that("uncapacitated model skips workload constraints", {
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

  dists <- surveyzones_compute_sparse_distances(pts, D_max = 10)

  # Uncapacitated: max_workload = Inf, target_zone_size = 3 -> K = 2
  plan <- surveyzones_build_zones(
    sparse_distances = dists,
    tracts = tracts,
    D_max = 10,
    target_zone_size = 3,
    enforce_partition = FALSE
  )

  expect_s3_class(plan, "surveyzones_plan")
  expect_equal(nrow(plan$zones), 2L)
  # No workload cap enforced — zones can have any workload
  expect_true(is.infinite(plan$parameters$max_workload_per_zone))
})

test_that("error when neither target_zone_size nor finite max_workload", {
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

  dists <- surveyzones_compute_sparse_distances(pts, D_max = 10)

  expect_error(
    surveyzones_build_zones(
      sparse_distances = dists,
      tracts = tracts,
      D_max = 10,
      enforce_partition = FALSE
    ),
    "target_zone_size"
  )
})

test_that("both target_zone_size and max_workload uses max K", {
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

  dists <- surveyzones_compute_sparse_distances(pts, D_max = 10)

  # target_zone_size = 3 -> K_size = 2
  # max_workload = 2 -> K_workload = ceil(6/2) = 3
  # Should use max(2, 3) = 3
  plan <- surveyzones_build_zones(
    sparse_distances = dists,
    tracts = tracts,
    D_max = 10,
    max_workload_per_zone = 2,
    target_zone_size = 3,
    enforce_partition = FALSE
  )

  expect_s3_class(plan, "surveyzones_plan")
  expect_gte(nrow(plan$zones), 3L)
  expect_true(all(plan$zones$total_workload <= 2))
})

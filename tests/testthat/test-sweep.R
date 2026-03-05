test_that("compare returns comparison tibble", {
  skip_if_not_installed("ROI.plugin.glpk")
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
    expected_service_time = rep(1, 4)
  )

  dists <- surveyzones_compute_sparse_distances(pts)

  plan_a <- surveyzones_build_zones(
    dists, tracts, D_max = 10, max_workload_per_zone = 2,
    enforce_partition = FALSE
  )
  plan_b <- surveyzones_build_zones(
    dists, tracts, D_max = 10, max_workload_per_zone = 4,
    enforce_partition = FALSE
  )

  result <- surveyzones_compare(list(tight = plan_a, loose = plan_b))

  expect_s3_class(result, "tbl_df")
  expect_equal(nrow(result), 2)
  expect_equal(result$plan_name, c("tight", "loose"))
})

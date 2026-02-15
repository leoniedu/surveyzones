test_that("sequencing returns a permutation of zone tracts", {
  # Create a small problem and solve it
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
    max_workload_per_zone = 4,
    enforce_partition = FALSE
  )

  # Sequence
  plan <- surveyzones_sequence(plan, dists)

  expect_false(is.null(plan$sequence))
  expect_true(all(c("zone_id", "tract_id", "visit_order") %in%
                    names(plan$sequence)))

  # Check each zone: sequenced tracts are a permutation of assigned tracts
  for (zid in unique(plan$zones$zone_id)) {
    assigned <- sort(plan$assignments$tract_id[
      plan$assignments$zone_id == zid
    ])
    sequenced <- sort(plan$sequence$tract_id[
      plan$sequence$zone_id == zid
    ])
    expect_equal(assigned, sequenced)
  }
})

test_that("sequence_zone handles 2-tract zones", {
  dists <- tibble::tibble(
    origin_id = c("A", "B"),
    destination_id = c("B", "A"),
    travel_time = c(5, 5)
  ) |> dplyr::arrange(origin_id, destination_id)

  result <- surveyzones_sequence_zone(c("A", "B"), dists)
  expect_equal(sort(result), c("A", "B"))
  expect_equal(length(result), 2)
})

test_that("method parameter is passed through", {
  dists <- tibble::tibble(
    origin_id = c("A", "B", "A", "C", "B", "C"),
    destination_id = c("B", "A", "C", "A", "C", "B"),
    travel_time = c(1, 1, 3, 3, 2, 2)
  ) |> dplyr::arrange(origin_id, destination_id)

  result <- surveyzones_sequence_zone(
    c("A", "B", "C"), dists, method = "nearest_insertion"
  )
  expect_equal(sort(result), c("A", "B", "C"))
  expect_equal(length(result), 3)
})

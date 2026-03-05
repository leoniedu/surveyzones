test_that("sequencing returns a permutation of zone tracts", {
  skip_if_not_installed("ROI.plugin.glpk")
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

  dists <- surveyzones_compute_sparse_distances(pts)

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
    distance = c(5, 5)
  ) |> dplyr::arrange(origin_id, destination_id)

  result <- surveyzones_sequence_zone(c("A", "B"), dists)
  expect_equal(sort(result), c("A", "B"))
  expect_equal(length(result), 2)
})

test_that("method parameter is passed through", {
  dists <- tibble::tibble(
    origin_id = c("A", "B", "A", "C", "B", "C"),
    destination_id = c("B", "A", "C", "A", "C", "B"),
    distance = c(1, 1, 3, 3, 2, 2)
  ) |> dplyr::arrange(origin_id, destination_id)

  result <- surveyzones_sequence_zone(
    c("A", "B", "C"), dists, method = "nearest_insertion"
  )
  expect_equal(sort(result), c("A", "B", "C"))
  expect_equal(length(result), 3)
})

test_that("MDS sequencing produces zone_score and zone_order", {
  # 4 zones in a line: Z1--Z2--Z3--Z4
  # Distances: adjacent = 1, skip-one = 2, skip-two = 3
  zones_tbl <- tibble::tibble(
    zone_id = paste0("Z", 1:4),
    partition_id = "P1",
    center_tract_id = paste0("C", 1:4),
    total_workload = 1,
    diameter = 0,
    n_tracts = 1L
  )

  # Complete symmetric distance matrix for 4 points on a line
  ids <- paste0("C", 1:4)
  pairs <- tidyr::expand_grid(origin_id = ids, destination_id = ids) |>
    dplyr::filter(origin_id != destination_id) |>
    dplyr::mutate(
      i = as.integer(gsub("C", "", origin_id)),
      j = as.integer(gsub("C", "", destination_id)),
      distance = abs(i - j)
    ) |>
    dplyr::select(origin_id, destination_id, distance)

  plan <- structure(
    list(
      zones = zones_tbl,
      assignments = tibble::tibble(
        tract_id = paste0("C", 1:4),
        zone_id = paste0("Z", 1:4),
        partition_id = "P1"
      ),
      zone_sequence = NULL,
      sequence = NULL,
      parameters = list(),
      diagnostics = list()
    ),
    class = "surveyzones_plan"
  )

  plan <- surveyzones_sequence_zones(plan, pairs, method = "mds")

  expect_false(is.null(plan$zone_sequence))
  expect_true("zone_score" %in% names(plan$zone_sequence))
  expect_true("zone_order" %in% names(plan$zone_sequence))
  expect_equal(nrow(plan$zone_sequence), 4)

  # zone_order should be 1:4
  expect_equal(sort(plan$zone_sequence$zone_order), 1:4)

  # MDS on a line should recover the linear order (up to reflection)
  scores <- plan$zone_sequence |>
    dplyr::arrange(zone_id) |>
    dplyr::pull(zone_score)

  # Consecutive scores should be monotonic (either all increasing or all decreasing)
  diffs <- diff(scores)
  expect_true(all(diffs > 0) || all(diffs < 0))
})

test_that("MDS sequencing handles single-zone partition", {
  zones_tbl <- tibble::tibble(
    zone_id = "Z1",
    partition_id = "P1",
    center_tract_id = "C1",
    total_workload = 1,
    diameter = 0,
    n_tracts = 1L
  )

  plan <- structure(
    list(
      zones = zones_tbl,
      assignments = tibble::tibble(
        tract_id = "C1", zone_id = "Z1", partition_id = "P1"
      ),
      zone_sequence = NULL,
      sequence = NULL,
      parameters = list(),
      diagnostics = list()
    ),
    class = "surveyzones_plan"
  )

  pairs <- tibble::tibble(
    origin_id = character(0),
    destination_id = character(0),
    distance = numeric(0)
  )

  plan <- surveyzones_sequence_zones(plan, pairs, method = "mds")
  expect_equal(plan$zone_sequence$zone_score, 0)
  expect_equal(plan$zone_sequence$zone_order, 1L)
})

test_that("MDS sequencing errors on incomplete distances", {
  zones_tbl <- tibble::tibble(
    zone_id = paste0("Z", 1:3),
    partition_id = "P1",
    center_tract_id = paste0("C", 1:3),
    total_workload = 1,
    diameter = 0,
    n_tracts = 1L
  )

  # Only C1-C2 pair, missing C1-C3 and C2-C3
  pairs <- tibble::tibble(
    origin_id = c("C1", "C2"),
    destination_id = c("C2", "C1"),
    distance = c(1, 1)
  )

  plan <- structure(
    list(
      zones = zones_tbl,
      assignments = tibble::tibble(
        tract_id = paste0("C", 1:3),
        zone_id = paste0("Z", 1:3),
        partition_id = "P1"
      ),
      zone_sequence = NULL,
      sequence = NULL,
      parameters = list(),
      diagnostics = list()
    ),
    class = "surveyzones_plan"
  )

  expect_error(
    surveyzones_sequence_zones(plan, pairs, method = "mds"),
    "missing"
  )
})

test_that("TSP sequencing returns NA zone_score", {
  zones_tbl <- tibble::tibble(
    zone_id = paste0("Z", 1:3),
    partition_id = "P1",
    center_tract_id = paste0("C", 1:3),
    total_workload = 1,
    diameter = 0,
    n_tracts = 1L
  )

  ids <- paste0("C", 1:3)
  pairs <- tidyr::expand_grid(origin_id = ids, destination_id = ids) |>
    dplyr::filter(origin_id != destination_id) |>
    dplyr::mutate(distance = 1) |>
    dplyr::select(origin_id, destination_id, distance)

  plan <- structure(
    list(
      zones = zones_tbl,
      assignments = tibble::tibble(
        tract_id = paste0("C", 1:3),
        zone_id = paste0("Z", 1:3),
        partition_id = "P1"
      ),
      zone_sequence = NULL,
      sequence = NULL,
      parameters = list(),
      diagnostics = list()
    ),
    class = "surveyzones_plan"
  )

  plan <- surveyzones_sequence_zones(plan, pairs, method = "nn")

  expect_true("zone_score" %in% names(plan$zone_sequence))
  expect_true(all(is.na(plan$zone_sequence$zone_score)))
})

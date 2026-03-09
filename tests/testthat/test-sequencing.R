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

test_that("sequence_tracts handles 2-tract zones", {
  dists <- tibble::tibble(
    origin_id = c("A", "B"),
    destination_id = c("B", "A"),
    distance = c(5, 5)
  ) |> dplyr::arrange(origin_id, destination_id)

  result <- surveyzones_sequence_tracts(c("A", "B"), dists)
  expect_equal(sort(result), c("A", "B"))
  expect_equal(length(result), 2)
})

test_that("TSP seriation method parameter is passed through", {
  dists <- tibble::tibble(
    origin_id = c("A", "B", "A", "C", "B", "C"),
    destination_id = c("B", "A", "C", "A", "C", "B"),
    distance = c(1, 1, 3, 3, 2, 2)
  ) |> dplyr::arrange(origin_id, destination_id)

  result <- surveyzones_sequence_tracts(
    c("A", "B", "C"), dists, method = "TSP"
  )
  expect_equal(sort(result), c("A", "B", "C"))
  expect_equal(length(result), 3)
})

test_that("zone sequencing produces correct zone_order for line topology", {
  # 4 zones in a line: Z1--Z2--Z3--Z4, all within threshold
  zones_tbl <- tibble::tibble(
    zone_id = paste0("Z", 1:4),
    partition_id = "P1",
    center_tract_id = paste0("C", 1:4),
    total_workload = 1,
    diameter = 0,
    n_tracts = 1L
  )

  ids <- paste0("C", 1:4)
  pairs <- tidyr::expand_grid(origin_id = ids, destination_id = ids) |>
    dplyr::filter(origin_id != destination_id) |>
    dplyr::mutate(
      i = as.integer(gsub("C", "", origin_id)),
      j = as.integer(gsub("C", "", destination_id)),
      distance = abs(i - j)
    ) |>
    dplyr::select(origin_id, destination_id, distance)

  pts <- sf::st_as_sf(
    data.frame(
      tract_id = ids,
      lon = -38.50 - (seq_along(ids) - 1) * 0.001,
      lat = -13.00 - (seq_along(ids) - 1) * 0.001
    ),
    coords = c("lon", "lat"), crs = 4326
  )

  plan <- structure(
    list(
      zones = zones_tbl,
      assignments = tibble::tibble(
        tract_id = paste0("C", 1:4),
        zone_id = paste0("Z", 1:4),
        partition_id = "P1"
      ),
      zone_sequence = NULL, sequence = NULL,
      parameters = list(), diagnostics = list()
    ),
    class = "surveyzones_plan"
  )

  plan <- surveyzones_sequence_zones(plan, pairs, pts, threshold = 100)

  expect_false(is.null(plan$zone_sequence))
  expect_true("zone_order" %in% names(plan$zone_sequence))
  expect_true("group_id" %in% names(plan$zone_sequence))
  expect_equal(nrow(plan$zone_sequence), 4)
  expect_equal(sort(plan$zone_sequence$zone_order), 1:4)

  # TSP on a line should recover the linear order (up to reflection)
  ordered_zones <- plan$zone_sequence |>
    dplyr::arrange(zone_order) |>
    dplyr::pull(zone_id)
  zone_nums <- as.integer(gsub("Z", "", ordered_zones))
  diffs <- diff(zone_nums)
  expect_true(all(diffs > 0) || all(diffs < 0))
})

test_that("seriation handles single-zone partition", {
  zones_tbl <- tibble::tibble(
    zone_id = "Z1",
    partition_id = "P1",
    center_tract_id = "C1",
    total_workload = 1,
    diameter = 0,
    n_tracts = 1L
  )

  pts <- sf::st_as_sf(
    data.frame(tract_id = "C1", lon = -38.50, lat = -13.00),
    coords = c("lon", "lat"), crs = 4326
  )

  plan <- structure(
    list(
      zones = zones_tbl,
      assignments = tibble::tibble(
        tract_id = "C1", zone_id = "Z1", partition_id = "P1"
      ),
      zone_sequence = NULL, sequence = NULL,
      parameters = list(), diagnostics = list()
    ),
    class = "surveyzones_plan"
  )

  pairs <- tibble::tibble(
    origin_id = character(0),
    destination_id = character(0),
    distance = numeric(0)
  )

  plan <- surveyzones_sequence_zones(plan, pairs, pts)
  expect_equal(plan$zone_sequence$zone_order, 1L)
})

test_that("Spectral tract sequencing produces correct visit_order", {
  # 5 tracts in a line within one zone
  ids <- paste0("T", 1:5)
  pairs <- tidyr::expand_grid(origin_id = ids, destination_id = ids) |>
    dplyr::filter(origin_id != destination_id) |>
    dplyr::mutate(
      i = as.integer(gsub("T", "", origin_id)),
      j = as.integer(gsub("T", "", destination_id)),
      distance = abs(i - j)
    ) |>
    dplyr::select(origin_id, destination_id, distance)

  plan <- structure(
    list(
      zones = tibble::tibble(
        zone_id = "Z1", partition_id = "P1",
        center_tract_id = "T1", total_workload = 5,
        diameter = 4, n_tracts = 5L
      ),
      assignments = tibble::tibble(
        tract_id = ids, zone_id = "Z1", partition_id = "P1"
      ),
      zone_sequence = NULL, sequence = NULL,
      parameters = list(), diagnostics = list()
    ),
    class = "surveyzones_plan"
  )

  plan <- surveyzones_sequence(plan, pairs, method = "Spectral")

  expect_false(is.null(plan$sequence))
  expect_true("visit_order" %in% names(plan$sequence))
  expect_equal(nrow(plan$sequence), 5)
  expect_equal(sort(plan$sequence$visit_order), 1:5)

  # Spectral on a line should recover linear order (up to reflection)
  ordered_tracts <- plan$sequence |>
    dplyr::arrange(visit_order) |>
    dplyr::pull(tract_id)
  tract_nums <- as.integer(gsub("T", "", ordered_tracts))
  diffs <- diff(tract_nums)
  expect_true(all(diffs > 0) || all(diffs < 0))
})

test_that("seriation handles single-tract zone", {
  plan <- structure(
    list(
      zones = tibble::tibble(
        zone_id = "Z1", partition_id = "P1",
        center_tract_id = "T1", total_workload = 1,
        diameter = 0, n_tracts = 1L
      ),
      assignments = tibble::tibble(
        tract_id = "T1", zone_id = "Z1", partition_id = "P1"
      ),
      zone_sequence = NULL, sequence = NULL,
      parameters = list(), diagnostics = list()
    ),
    class = "surveyzones_plan"
  )

  pairs <- tibble::tibble(
    origin_id = character(0),
    destination_id = character(0),
    distance = numeric(0)
  )

  plan <- surveyzones_sequence(plan, pairs, method = "Spectral")
  expect_equal(plan$sequence$visit_order, 1L)
})

test_that("zone sequencing completes incomplete distances via haversine", {
  # 3 zone centers, but only provide C1-C2 distances
  zones_tbl <- tibble::tibble(
    zone_id = paste0("Z", 1:3),
    partition_id = "P1",
    center_tract_id = paste0("C", 1:3),
    total_workload = 1,
    diameter = 0,
    n_tracts = 1L
  )

  incomplete_pairs <- tibble::tibble(
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
      zone_sequence = NULL, sequence = NULL,
      parameters = list(), diagnostics = list()
    ),
    class = "surveyzones_plan"
  )

  pts <- sf::st_as_sf(
    data.frame(
      tract_id = paste0("C", 1:3),
      lon = c(-38.50, -38.501, -38.502),
      lat = c(-13.00, -13.001, -13.002)
    ),
    coords = c("lon", "lat"), crs = 4326
  )

  # Should succeed — missing pairs filled via haversine
  plan <- surveyzones_sequence_zones(plan, incomplete_pairs, pts)
  expect_false(is.null(plan$zone_sequence))
  expect_equal(nrow(plan$zone_sequence), 3)
  expect_equal(sort(plan$zone_sequence$zone_order), 1:3)
})

test_that("by_partition = FALSE sequences all zones together", {
  zones_tbl <- tibble::tibble(
    zone_id = paste0("Z", 1:4),
    partition_id = c("P1", "P1", "P2", "P2"),
    center_tract_id = paste0("C", 1:4),
    total_workload = 1,
    diameter = 0,
    n_tracts = 1L
  )

  ids <- paste0("C", 1:4)
  pairs <- tidyr::expand_grid(origin_id = ids, destination_id = ids) |>
    dplyr::filter(origin_id != destination_id) |>
    dplyr::mutate(
      i = as.integer(gsub("C", "", origin_id)),
      j = as.integer(gsub("C", "", destination_id)),
      distance = abs(i - j)
    ) |>
    dplyr::select(origin_id, destination_id, distance)

  pts <- sf::st_as_sf(
    data.frame(
      tract_id = ids,
      lon = -38.50 - (seq_along(ids) - 1) * 0.001,
      lat = -13.00 - (seq_along(ids) - 1) * 0.001
    ),
    coords = c("lon", "lat"), crs = 4326
  )

  plan <- structure(
    list(
      zones = zones_tbl,
      assignments = tibble::tibble(
        tract_id = paste0("C", 1:4),
        zone_id = paste0("Z", 1:4),
        partition_id = c("P1", "P1", "P2", "P2")
      ),
      zone_sequence = NULL, sequence = NULL,
      parameters = list(), diagnostics = list()
    ),
    class = "surveyzones_plan"
  )

  plan <- surveyzones_sequence_zones(
    plan, pairs, pts, threshold = 100, by_partition = FALSE
  )

  expect_false(is.null(plan$zone_sequence))
  expect_equal(nrow(plan$zone_sequence), 4)
  # All 4 zones get a single global ordering 1:4
  expect_equal(sort(plan$zone_sequence$zone_order), 1:4)
  # partition_id column still present
  expect_true("partition_id" %in% names(plan$zone_sequence))
})

test_that("OLO seriation method works", {
  ids <- paste0("T", 1:4)
  pairs <- tidyr::expand_grid(origin_id = ids, destination_id = ids) |>
    dplyr::filter(origin_id != destination_id) |>
    dplyr::mutate(
      i = as.integer(gsub("T", "", origin_id)),
      j = as.integer(gsub("T", "", destination_id)),
      distance = abs(i - j)
    ) |>
    dplyr::select(origin_id, destination_id, distance)

  result <- surveyzones_sequence_tracts(ids, pairs, method = "OLO")
  expect_equal(sort(result), sort(ids))
  expect_equal(length(result), 4)
})

test_that(".build_zone_dist_matrix computes min pairwise distances", {
  # 4 tracts in 2 zones: Z1 = {T1, T2}, Z2 = {T3, T4}
  # T2 and T3 are close (distance 1), T1-T4 are far
  assignments <- tibble::tibble(
    tract_id = paste0("T", 1:4),
    zone_id = c("Z1", "Z1", "Z2", "Z2"),
    partition_id = "P1"
  )

  sparse_dists <- tibble::tibble(
    origin_id = c("T1", "T2", "T1", "T3", "T2", "T4", "T3", "T4",
                  "T2", "T1", "T3", "T1", "T4", "T2", "T4", "T3"),
    destination_id = c("T2", "T1", "T3", "T1", "T3", "T2", "T4", "T3",
                       "T4", "T4", "T2", "T2", "T1", "T3", "T1", "T4"),
    distance = c(2, 2, 10, 10, 1, 8, 3, 3,
                 8, 9, 1, 2, 9, 3, 10, 3)
  )
  # Cross-zone pairs: T1-T3=10, T2-T3=1, T1-T4=9, T2-T4=8
  # min(Z1, Z2) should be 1 (from T2-T3)

  mat <- surveyzones:::.build_zone_dist_matrix(
    zone_ids = c("Z1", "Z2"),
    assignments = assignments,
    sparse_distances = sparse_dists
  )

  expect_equal(nrow(mat), 2)
  expect_equal(ncol(mat), 2)
  expect_equal(rownames(mat), c("Z1", "Z2"))
  expect_equal(mat["Z1", "Z2"], 1)
  expect_equal(mat["Z2", "Z1"], 1)
  expect_equal(mat["Z1", "Z1"], 0)
})

test_that(".build_zone_dist_matrix has no Inf when distances are complete", {
  # 3 zones with 2 tracts each, complete distance table
  assignments <- tibble::tibble(
    tract_id = paste0("T", 1:6),
    zone_id = rep(paste0("Z", 1:3), each = 2),
    partition_id = "P1"
  )

  ids <- paste0("T", 1:6)
  pairs <- tidyr::expand_grid(origin_id = ids, destination_id = ids) |>
    dplyr::filter(origin_id != destination_id) |>
    dplyr::mutate(
      i = as.integer(gsub("T", "", origin_id)),
      j = as.integer(gsub("T", "", destination_id)),
      distance = abs(i - j)
    ) |>
    dplyr::select(origin_id, destination_id, distance)

  mat <- surveyzones:::.build_zone_dist_matrix(
    zone_ids = paste0("Z", 1:3),
    assignments = assignments,
    sparse_distances = pairs
  )

  # No Inf entries when all pairs are present
  expect_true(all(is.finite(mat)))
  # Diagonal is 0
  expect_equal(diag(mat), c(Z1 = 0, Z2 = 0, Z3 = 0))
  # Symmetric
  expect_equal(mat, t(mat))
})

test_that("sequence_zones uses min pairwise distances, not center-to-center", {
  # Z1 center at C1 (far left), Z2 center at C4 (far right)
  # But Z1 contains T2 (near right) and Z2 contains T3 (near left)
  # So min-pairwise(Z1,Z2) = d(T2,T3) which is small
  # Z3 center at C5 (middle), but all its tracts are in the middle
  #
  # Layout: C1...T2|T3...C4   C5
  # Center-to-center: d(Z1,Z3) < d(Z1,Z2) because C5 is closer to C1 than C4
  # Min-pairwise: d(Z1,Z2) < d(Z1,Z3) because T2-T3 are touching

  assignments <- tibble::tibble(
    tract_id = c("C1", "T2", "T3", "C4", "C5", "T6"),
    zone_id = c("Z1", "Z1", "Z2", "Z2", "Z3", "Z3"),
    partition_id = "P1"
  )

  zones_tbl <- tibble::tibble(
    zone_id = c("Z1", "Z2", "Z3"),
    partition_id = "P1",
    center_tract_id = c("C1", "C4", "C5"),
    total_workload = 2,
    diameter = 1,
    n_tracts = 2L
  )

  # All 6 tracts along a line: C1(0), T2(4), T3(5), C4(9), C5(15), T6(16)
  positions <- c(C1 = 0, T2 = 4, T3 = 5, C4 = 9, C5 = 15, T6 = 16)
  ids <- names(positions)
  pairs <- tidyr::expand_grid(origin_id = ids, destination_id = ids) |>
    dplyr::filter(origin_id != destination_id) |>
    dplyr::mutate(
      distance = abs(positions[origin_id] - positions[destination_id])
    )

  pts <- sf::st_as_sf(
    data.frame(
      tract_id = ids,
      lon = -38.50 - (positions / 1000),
      lat = rep(-13.00, length(ids))
    ),
    coords = c("lon", "lat"), crs = 4326
  )

  plan <- structure(
    list(
      zones = zones_tbl,
      assignments = assignments,
      zone_sequence = NULL, sequence = NULL,
      parameters = list(), diagnostics = list(),
      access_points = pts
    ),
    class = "surveyzones_plan"
  )

  plan <- surveyzones_sequence_zones(
    plan, pairs, pts,
    threshold = 100, complete_distances = FALSE
  )

  # Z1 and Z2 should be adjacent in the ordering (min pairwise = 1)
  # Z3 should be at one end (far from both)
  ordered <- plan$zone_sequence |>
    dplyr::arrange(zone_order)

  z1_order <- ordered$zone_order[ordered$zone_id == "Z1"]
  z2_order <- ordered$zone_order[ordered$zone_id == "Z2"]

  # Z1 and Z2 must be consecutive
  expect_equal(abs(z1_order - z2_order), 1)
})

test_that("zone sequencing groups zones by threshold (gap-based)", {
  # 6 zones in two clusters: Z1-Z3 close together, Z4-Z6 close together
  # Clusters are far apart (distance 50)
  zones_tbl <- tibble::tibble(
    zone_id = paste0("Z", 1:6),
    partition_id = "P1",
    center_tract_id = paste0("C", 1:6),
    total_workload = 1,
    diameter = 0,
    n_tracts = 1L
  )

  # C1(0), C2(2), C3(4), C4(54), C5(56), C6(58)
  positions <- c(C1 = 0, C2 = 2, C3 = 4, C4 = 54, C5 = 56, C6 = 58)
  ids <- names(positions)
  pairs <- tidyr::expand_grid(origin_id = ids, destination_id = ids) |>
    dplyr::filter(origin_id != destination_id) |>
    dplyr::mutate(
      distance = abs(positions[origin_id] - positions[destination_id])
    )

  pts <- sf::st_as_sf(
    data.frame(
      tract_id = ids,
      lon = -38.50 - (positions / 10000),
      lat = rep(-13.00, length(ids))
    ),
    coords = c("lon", "lat"), crs = 4326
  )

  plan <- structure(
    list(
      zones = zones_tbl,
      assignments = tibble::tibble(
        tract_id = ids,
        zone_id = paste0("Z", 1:6),
        partition_id = "P1"
      ),
      zone_sequence = NULL, sequence = NULL,
      parameters = list(), diagnostics = list()
    ),
    class = "surveyzones_plan"
  )

  plan <- surveyzones_sequence_zones(
    plan, pairs, pts, threshold = 10, complete_distances = FALSE
  )

  ordered <- plan$zone_sequence |> dplyr::arrange(zone_order)

  # Z1-Z3 should be consecutive and Z4-Z6 should be consecutive
  cluster1_orders <- ordered$zone_order[ordered$zone_id %in% paste0("Z", 1:3)]
  cluster2_orders <- ordered$zone_order[ordered$zone_id %in% paste0("Z", 4:6)]

  # Each cluster occupies 3 consecutive positions
  expect_equal(diff(sort(cluster1_orders)), c(1L, 1L))
  expect_equal(diff(sort(cluster2_orders)), c(1L, 1L))

  # The two clusters don't interleave
  expect_true(
    max(cluster1_orders) < min(cluster2_orders) ||
    max(cluster2_orders) < min(cluster1_orders)
  )

  # group_id should distinguish the two clusters
  expect_true("group_id" %in% names(ordered))
  cluster1_groups <- unique(ordered$group_id[ordered$zone_id %in% paste0("Z", 1:3)])
  cluster2_groups <- unique(ordered$group_id[ordered$zone_id %in% paste0("Z", 4:6)])
  expect_length(cluster1_groups, 1)
  expect_length(cluster2_groups, 1)
  expect_false(cluster1_groups == cluster2_groups)
})

test_that(".orient_tract_sequences reverses tracts when exit is closer to last", {
  # Zone A has tracts T1 -> T2 -> T3 (visit_order 1,2,3)
  # Zone B has tracts T4 -> T5 -> T6 (visit_order 1,2,3)
  # Exit of zone A = T3, entry of zone B = T4, end of zone B = T6
  # If T3 is closer to T6 than T4, zone B should be reversed
  sequence <- tibble::tibble(
    zone_id  = c(rep("ZA", 3), rep("ZB", 3)),
    tract_id = c("T1", "T2", "T3", "T4", "T5", "T6"),
    visit_order = c(1L, 2L, 3L, 1L, 2L, 3L)
  )

  zone_sequence <- tibble::tibble(
    partition_id = "P1",
    zone_id = c("ZA", "ZB"),
    zone_order = 1:2,
    group_id = c(1L, 1L)
  )

  # T3 -> T4 = 10, T3 -> T6 = 2 (closer to last, should reverse)
  sparse_distances <- tibble::tibble(
    origin_id      = c("T3", "T3"),
    destination_id = c("T4", "T6"),
    distance       = c(10,   2)
  )

  result <- surveyzones:::.orient_tract_sequences(
    sequence, zone_sequence, sparse_distances
  )

  zb <- result |> dplyr::filter(zone_id == "ZB") |> dplyr::arrange(visit_order)
  # After reversal: T6 should be first (visit_order 1), T4 last (visit_order 3)
  expect_equal(zb$tract_id, c("T6", "T5", "T4"))

  # Zone A should be unchanged
  za <- result |> dplyr::filter(zone_id == "ZA") |> dplyr::arrange(visit_order)
  expect_equal(za$tract_id, c("T1", "T2", "T3"))
})

test_that(".orient_tract_sequences does not reverse when entry is already closer", {
  sequence <- tibble::tibble(
    zone_id  = c(rep("ZA", 3), rep("ZB", 3)),
    tract_id = c("T1", "T2", "T3", "T4", "T5", "T6"),
    visit_order = c(1L, 2L, 3L, 1L, 2L, 3L)
  )

  zone_sequence <- tibble::tibble(
    partition_id = "P1",
    zone_id = c("ZA", "ZB"),
    zone_order = 1:2,
    group_id = c(1L, 1L)
  )

  # T3 -> T4 = 2 (closer to first, no reversal), T3 -> T6 = 10
  sparse_distances <- tibble::tibble(
    origin_id      = c("T3", "T3"),
    destination_id = c("T4", "T6"),
    distance       = c(2,   10)
  )

  result <- surveyzones:::.orient_tract_sequences(
    sequence, zone_sequence, sparse_distances
  )

  zb <- result |> dplyr::filter(zone_id == "ZB") |> dplyr::arrange(visit_order)
  expect_equal(zb$tract_id, c("T4", "T5", "T6"))
})

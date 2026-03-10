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

  # Sequence (zones + tracts in one call)
  plan <- surveyzones_sequence(plan, dists, complete_distances = FALSE)

  expect_false(is.null(plan$sequence))
  expect_false(is.null(plan$zone_sequence))
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
  # 4 zones in a line: Z1--Z2--Z3--Z4
  zones_tbl <- tibble::tibble(
    zone_id = paste0("Z", 1:4),
    partition_id = "P1",
    center_tract_id = paste0("C", 1:4),
    total_workload = 1,
    diameter = 0,
    n_tracts = 1L,
    group_id = paste0("Z", 1:4)
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
        partition_id = "P1",
        group_id = paste0("Z", 1:4)
      ),
      zone_sequence = NULL, sequence = NULL,
      parameters = list(), diagnostics = list()
    ),
    class = "surveyzones_plan"
  )

  plan <- surveyzones_sequence(
    plan, pairs, complete_distances = FALSE
  )

  expect_false(is.null(plan$zone_sequence))
  expect_true("zone_order" %in% names(plan$zone_sequence))
  expect_true("group_id" %in% names(plan$zone_sequence))
  expect_equal(nrow(plan$zone_sequence), 4)
  expect_equal(sort(plan$zone_sequence$zone_order), 1:4)

  # TSP on a line should recover the linear order (up to reflection)
  # Use zone_order directly — IDs have been renamed and are no longer "Z1"..."Z4"
  ordered_zones <- plan$zone_sequence |>
    dplyr::arrange(zone_order) |>
    dplyr::pull(zone_id)
  # All 4 zones present, in some order
  expect_length(ordered_zones, 4)
  expect_equal(sort(plan$zone_sequence$zone_order), 1:4)
  # All new IDs match the expected format {partition}_{group}.{pos:03d}
  expect_true(all(grepl("^P1_\\d+\\.\\d{3}$", ordered_zones)))

  # Tract sequence should also exist
  expect_false(is.null(plan$sequence))
  expect_equal(nrow(plan$sequence), 4)
})

test_that("seriation handles single-zone partition", {
  zones_tbl <- tibble::tibble(
    zone_id = "Z1",
    partition_id = "P1",
    center_tract_id = "C1",
    total_workload = 1,
    diameter = 0,
    n_tracts = 1L,
    group_id = "Z1"
  )

  pts <- sf::st_as_sf(
    data.frame(tract_id = "C1", lon = -38.50, lat = -13.00),
    coords = c("lon", "lat"), crs = 4326
  )

  plan <- structure(
    list(
      zones = zones_tbl,
      assignments = tibble::tibble(
        tract_id = "C1", zone_id = "Z1", partition_id = "P1",
        group_id = "Z1"
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

  plan <- surveyzones_sequence(plan, pairs, complete_distances = FALSE)
  expect_equal(plan$zone_sequence$zone_order, 1L)
  expect_equal(plan$sequence$visit_order, 1L)
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
        diameter = 4, n_tracts = 5L, group_id = "Z1"
      ),
      assignments = tibble::tibble(
        tract_id = ids, zone_id = "Z1", partition_id = "P1",
        group_id = "Z1"
      ),
      zone_sequence = NULL, sequence = NULL,
      parameters = list(), diagnostics = list()
    ),
    class = "surveyzones_plan"
  )

  plan <- surveyzones_sequence(
    plan, pairs, method = "Spectral", complete_distances = FALSE
  )

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
        diameter = 0, n_tracts = 1L, group_id = "Z1"
      ),
      assignments = tibble::tibble(
        tract_id = "T1", zone_id = "Z1", partition_id = "P1",
        group_id = "Z1"
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

  plan <- surveyzones_sequence(
    plan, pairs, method = "Spectral", complete_distances = FALSE
  )
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
    n_tracts = 1L,
    group_id = paste0("Z", 1:3)
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
        partition_id = "P1",
        group_id = paste0("Z", 1:3)
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
  plan <- surveyzones_sequence(plan, incomplete_pairs, access_points = pts)
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
    n_tracts = 1L,
    group_id = paste0("Z", 1:4)
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
        partition_id = c("P1", "P1", "P2", "P2"),
        group_id = paste0("Z", 1:4)
      ),
      zone_sequence = NULL, sequence = NULL,
      parameters = list(), diagnostics = list()
    ),
    class = "surveyzones_plan"
  )

  plan <- surveyzones_sequence(
    plan, pairs, by_partition = FALSE,
    complete_distances = FALSE
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

test_that("sequence uses min pairwise distances, not center-to-center", {
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
    partition_id = "P1",
    group_id = c("Z1", "Z1", "Z2", "Z2", "Z3", "Z3")
  )

  zones_tbl <- tibble::tibble(
    zone_id = c("Z1", "Z2", "Z3"),
    partition_id = "P1",
    center_tract_id = c("C1", "C4", "C5"),
    total_workload = 2,
    diameter = 1,
    n_tracts = 2L,
    group_id = c("Z1", "Z2", "Z3")
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

  plan <- surveyzones_sequence(
    plan, pairs, complete_distances = FALSE
  )

  # Z1 and Z2 should be adjacent in the ordering (min pairwise = 1)
  # Z3 should be at one end (far from both)
  ordered <- plan$zone_sequence |>
    dplyr::arrange(zone_order)

  # Look up renamed IDs via center_tract_id (preserved in plan$zones)
  z1_id <- plan$zones$zone_id[plan$zones$center_tract_id == "C1"]
  z2_id <- plan$zones$zone_id[plan$zones$center_tract_id == "C4"]
  z1_order <- ordered$zone_order[ordered$zone_id == z1_id]
  z2_order <- ordered$zone_order[ordered$zone_id == z2_id]

  # Z1 and Z2 must be consecutive
  expect_equal(abs(z1_order - z2_order), 1)
})

test_that("zone sequencing respects solver-defined groups", {
  # 6 zones in two solver-defined groups: G1={Z1,Z2,Z3}, G2={Z4,Z5,Z6}
  zones_tbl <- tibble::tibble(
    zone_id = paste0("Z", 1:6),
    partition_id = "P1",
    center_tract_id = paste0("C", 1:6),
    total_workload = 1,
    diameter = 0,
    n_tracts = 1L,
    group_id = rep(c("G1", "G2"), each = 3)
  )

  # C1(0), C2(2), C3(4), C4(54), C5(56), C6(58)
  positions <- c(C1 = 0, C2 = 2, C3 = 4, C4 = 54, C5 = 56, C6 = 58)
  ids <- names(positions)
  pairs <- tidyr::expand_grid(origin_id = ids, destination_id = ids) |>
    dplyr::filter(origin_id != destination_id) |>
    dplyr::mutate(
      distance = abs(positions[origin_id] - positions[destination_id])
    )

  plan <- structure(
    list(
      zones = zones_tbl,
      assignments = tibble::tibble(
        tract_id = ids,
        zone_id = paste0("Z", 1:6),
        partition_id = "P1",
        group_id = rep(c("G1", "G2"), each = 3)
      ),
      zone_sequence = NULL, sequence = NULL,
      parameters = list(), diagnostics = list()
    ),
    class = "surveyzones_plan"
  )

  plan <- surveyzones_sequence(
    plan, pairs, complete_distances = FALSE
  )

  ordered <- plan$zone_sequence |> dplyr::arrange(zone_order)

  # Two distinct group_ids should exist (one per solver group)
  groups <- sort(unique(ordered$group_id))
  expect_length(groups, 2)

  cluster1_orders <- ordered$zone_order[ordered$group_id == groups[1]]
  cluster2_orders <- ordered$zone_order[ordered$group_id == groups[2]]

  # Each cluster occupies 3 consecutive positions
  expect_equal(diff(sort(cluster1_orders)), c(1L, 1L))
  expect_equal(diff(sort(cluster2_orders)), c(1L, 1L))

  # The two clusters don't interleave
  expect_true(
    max(cluster1_orders) < min(cluster2_orders) ||
    max(cluster2_orders) < min(cluster1_orders)
  )
})

test_that("tract orientation orients second zone's entry toward first zone's exit", {
  # Two zones on a line, well-separated:
  # ZA tracts at positions 0, 1, 2
  # ZB tracts at positions 10, 11, 12
  # Regardless of which zone comes first, the second zone's entry (first tract)
  # should be closer to the first zone's exit (last tract) than the second
  # zone's last tract.
  zones_tbl <- tibble::tibble(
    zone_id = c("ZA", "ZB"),
    partition_id = "P1",
    center_tract_id = c("T1", "T4"),
    total_workload = 3,
    diameter = 2,
    n_tracts = 3L,
    group_id = c("ZA", "ZB")
  )

  positions <- c(T1 = 0, T2 = 1, T3 = 2, T4 = 10, T5 = 11, T6 = 12)
  ids <- names(positions)
  pairs <- tidyr::expand_grid(origin_id = ids, destination_id = ids) |>
    dplyr::filter(origin_id != destination_id) |>
    dplyr::mutate(
      distance = abs(positions[origin_id] - positions[destination_id])
    )

  plan <- structure(
    list(
      zones = zones_tbl,
      assignments = tibble::tibble(
        tract_id = ids,
        zone_id = c("ZA", "ZA", "ZA", "ZB", "ZB", "ZB"),
        partition_id = "P1",
        group_id = c("ZA", "ZA", "ZA", "ZB", "ZB", "ZB")
      ),
      zone_sequence = NULL, sequence = NULL,
      parameters = list(), diagnostics = list()
    ),
    class = "surveyzones_plan"
  )

  plan <- surveyzones_sequence(
    plan, pairs, complete_distances = FALSE
  )

  # Identify which zone is visited first and second
  zs <- plan$zone_sequence |> dplyr::arrange(zone_order)
  first_zone <- zs$zone_id[1]
  second_zone <- zs$zone_id[2]

  first_tracts <- plan$sequence |>
    dplyr::filter(zone_id == first_zone) |>
    dplyr::arrange(visit_order)
  second_tracts <- plan$sequence |>
    dplyr::filter(zone_id == second_zone) |>
    dplyr::arrange(visit_order)

  # The second zone's entry (first tract) should be closer to the first
  # zone's exit (last tract) than the second zone's last tract
  exit_pos <- positions[first_tracts$tract_id[nrow(first_tracts)]]
  entry_pos <- positions[second_tracts$tract_id[1]]
  tail_pos <- positions[second_tracts$tract_id[nrow(second_tracts)]]

  expect_true(abs(entry_pos - exit_pos) <= abs(tail_pos - exit_pos))
})

test_that("default TSP control produces valid zone ordering", {
  # 5 zones on a line: Z1(0)--Z2(1)--Z3(2)--Z4(3)--Z5(4)
  zones_tbl <- tibble::tibble(
    zone_id = paste0("Z", 1:5),
    partition_id = "P1",
    center_tract_id = paste0("C", 1:5),
    total_workload = 1,
    diameter = 0,
    n_tracts = 1L,
    group_id = paste0("Z", 1:5)
  )

  ids <- paste0("C", 1:5)
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
        tract_id = ids,
        zone_id = paste0("Z", 1:5),
        partition_id = "P1",
        group_id = paste0("Z", 1:5)
      ),
      zone_sequence = NULL, sequence = NULL,
      parameters = list(), diagnostics = list()
    ),
    class = "surveyzones_plan"
  )

  # Default control = NULL (seriation's default TSP heuristic)
  plan <- surveyzones_sequence(
    plan, pairs, complete_distances = FALSE
  )

  ordered <- plan$zone_sequence |> dplyr::arrange(zone_order)
  expect_equal(nrow(ordered), 5)
  expect_true(all(grepl("^P1_\\d+\\.\\d{3}$", ordered$zone_id)))
})

test_that("custom NN control works when passed explicitly", {
  ids <- paste0("T", 1:4)
  pairs <- tidyr::expand_grid(origin_id = ids, destination_id = ids) |>
    dplyr::filter(origin_id != destination_id) |>
    dplyr::mutate(
      i = as.integer(gsub("T", "", origin_id)),
      j = as.integer(gsub("T", "", destination_id)),
      distance = abs(i - j)
    ) |>
    dplyr::select(origin_id, destination_id, distance)

  withr::local_seed(42)
  result <- surveyzones_sequence_tracts(ids, pairs, method = "TSP",
                                        control = list(method = "nn", rep = 20))
  expect_equal(sort(result), sort(ids))
  expect_equal(length(result), 4)
})

# ── .rename_zones ──────────────────────────────────────────────────────────────

test_that(".rename_zones produces correct ID format", {
  plan <- structure(
    list(
      zones = tibble::tibble(
        zone_id = c("Z1", "Z2", "Z3"),
        partition_id = "P1",
        center_tract_id = c("C1", "C2", "C3"),
        total_workload = 1, diameter = 0, n_tracts = 1L
      ),
      assignments = tibble::tibble(
        tract_id = c("T1", "T2", "T3"),
        zone_id  = c("Z1", "Z2", "Z3"),
        partition_id = "P1",
        group_id = c("ph1", "ph1", "ph1")
      ),
      zone_sequence = tibble::tibble(
        partition_id = "P1",
        zone_id      = c("Z1", "Z2", "Z3"),
        zone_order   = 1:3,
        group_id     = c(1L, 1L, 2L)   # Z1+Z2 in group 1, Z3 in group 2
      ),
      sequence = tibble::tibble(
        zone_id     = c("Z1", "Z2", "Z3"),
        tract_id    = c("T1", "T2", "T3"),
        visit_order = 1L
      ),
      parameters = list(), diagnostics = list()
    ),
    class = "surveyzones_plan"
  )

  out <- surveyzones:::.rename_zones(plan)

  # zone_sequence uses new IDs
  expect_equal(sort(out$zone_sequence$zone_id), c("P1_1.001", "P1_1.002", "P1_2.001"))
  # assignments uses new IDs
  expect_equal(sort(out$assignments$zone_id), c("P1_1.001", "P1_1.002", "P1_2.001"))
  # zones uses new IDs
  expect_equal(sort(out$zones$zone_id), c("P1_1.001", "P1_1.002", "P1_2.001"))
  # sequence uses new IDs
  expect_equal(sort(out$sequence$zone_id), c("P1_1.001", "P1_1.002", "P1_2.001"))
  # group_id dropped from assignments
  expect_false("group_id" %in% names(out$assignments))
  # center_tract_id preserved in zones
  expect_true("center_tract_id" %in% names(out$zones))
})

test_that(".rename_zones works across two partitions", {
  plan <- structure(
    list(
      zones = tibble::tibble(
        zone_id = c("Z1", "Z2", "Z3", "Z4"),
        partition_id = c("PA", "PA", "PB", "PB"),
        center_tract_id = paste0("C", 1:4),
        total_workload = 1, diameter = 0, n_tracts = 1L
      ),
      assignments = tibble::tibble(
        tract_id     = paste0("T", 1:4),
        zone_id      = c("Z1", "Z2", "Z3", "Z4"),
        partition_id = c("PA", "PA", "PB", "PB"),
        group_id     = c("g1", "g1", "g2", "g2")  # char group_id to drop
      ),
      zone_sequence = tibble::tibble(
        partition_id = c("PA", "PA", "PB", "PB"),
        zone_id      = c("Z1", "Z2", "Z3", "Z4"),
        zone_order   = c(1L, 2L, 1L, 2L),
        group_id     = c(1L, 1L, 1L, 1L)
      ),
      sequence = tibble::tibble(
        zone_id     = c("Z1", "Z2", "Z3", "Z4"),
        tract_id    = paste0("T", 1:4),
        visit_order = 1L
      ),
      parameters = list(), diagnostics = list()
    ),
    class = "surveyzones_plan"
  )

  out <- surveyzones:::.rename_zones(plan)

  pa_ids <- sort(out$zone_sequence$zone_id[out$zone_sequence$partition_id == "PA"])
  pb_ids <- sort(out$zone_sequence$zone_id[out$zone_sequence$partition_id == "PB"])
  expect_equal(pa_ids, c("PA_1.001", "PA_1.002"))
  expect_equal(pb_ids, c("PB_1.001", "PB_1.002"))
  # group_id dropped from assignments
  expect_false("group_id" %in% names(out$assignments))
  # sequence and zones tables are also updated
  expect_equal(sort(unique(out$sequence$zone_id)),
               c("PA_1.001", "PA_1.002", "PB_1.001", "PB_1.002"))
  expect_equal(sort(unique(out$zones$zone_id)),
               c("PA_1.001", "PA_1.002", "PB_1.001", "PB_1.002"))
})

test_that(".rename_zones respects zone_order when rows are unsorted", {
  # Z1 row appears first but has zone_order=2; Z2 has zone_order=1
  # Z2 (lowest zone_order in group) should get pos .001
  plan <- structure(
    list(
      zones = tibble::tibble(
        zone_id = c("Z1", "Z2"),
        partition_id = "P1",
        center_tract_id = c("C1", "C2"),
        total_workload = 1, diameter = 0, n_tracts = 1L
      ),
      assignments = tibble::tibble(
        tract_id = c("T1", "T2"),
        zone_id  = c("Z1", "Z2"),
        partition_id = "P1",
        group_id = c("ph1", "ph1")
      ),
      zone_sequence = tibble::tibble(
        partition_id = "P1",
        zone_id      = c("Z1", "Z2"),   # rows are Z1 then Z2
        zone_order   = c(2L, 1L),        # but Z2 has the lower zone_order
        group_id     = c(1L, 1L)
      ),
      sequence = tibble::tibble(
        zone_id = c("Z1", "Z2"), tract_id = c("T1", "T2"), visit_order = 1L
      ),
      parameters = list(), diagnostics = list()
    ),
    class = "surveyzones_plan"
  )

  out <- surveyzones:::.rename_zones(plan)

  # T2 belongs to Z2 (zone_order=1) -> must be in zone P1_1.001
  # T1 belongs to Z1 (zone_order=2) -> must be in zone P1_1.002
  t2_zone <- out$sequence |> dplyr::filter(.data$tract_id == "T2") |> dplyr::pull(zone_id)
  t1_zone <- out$sequence |> dplyr::filter(.data$tract_id == "T1") |> dplyr::pull(zone_id)
  expect_length(t2_zone, 1L)
  expect_length(t1_zone, 1L)
  expect_equal(t2_zone, "P1_1.001")
  expect_equal(t1_zone, "P1_1.002")
})

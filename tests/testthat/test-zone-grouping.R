test_that(".validate_zone_assignments catches missing columns", {
  expect_error(
    surveyzones:::.validate_zone_assignments(data.frame(x = 1)),
    "missing column"
  )
})

test_that(".validate_zone_assignments catches duplicate tract_id across zones", {
  za <- data.frame(
    tract_id = c("T1", "T1", "T2"),
    zone_id  = c("Z1", "Z2", "Z2")
  )
  expect_error(
    surveyzones:::.validate_zone_assignments(za),
    "must be unique"
  )
})

test_that(".validate_zone_assignments catches zones spanning partitions", {
  za <- data.frame(
    tract_id = c("T1", "T2", "T3"),
    zone_id  = c("Z1", "Z1", "Z2"),
    partition_id = c("A", "B", "B")
  )
  expect_error(
    surveyzones:::.validate_zone_assignments(za),
    "span more than one partition"
  )
})

test_that(".validate_zone_assignments defaults partition_id to 'all'", {
  za <- data.frame(tract_id = c("T1", "T2"), zone_id = c("Z1", "Z1"))
  out <- surveyzones:::.validate_zone_assignments(za)
  expect_equal(unique(out$partition_id), "all")
})

test_that(".validate_zone_assignments passes through existing partition_id", {
  za <- data.frame(
    tract_id = c("T1", "T2"), zone_id = c("Z1", "Z2"),
    partition_id = c("A", "B")
  )
  out <- surveyzones:::.validate_zone_assignments(za)
  expect_equal(out$partition_id, c("A", "B"))
})

test_that(".find_zone_centers picks the medoid minimizing total distance", {
  za <- surveyzones:::.validate_zone_assignments(data.frame(
    tract_id = c("T1", "T2", "T3"),
    zone_id  = c("Z1", "Z1", "Z1")
  ))
  # Positions: T1=0, T2=1, T3=3 -> sums: T1=4, T2=3, T3=5 -> medoid T2
  positions <- c(T1 = 0, T2 = 1, T3 = 3)
  ids <- names(positions)
  dists <- tidyr::expand_grid(origin_id = ids, destination_id = ids) |>
    dplyr::filter(origin_id != destination_id) |>
    dplyr::mutate(distance = abs(positions[origin_id] - positions[destination_id]))

  centers <- surveyzones:::.find_zone_centers(za, dists)
  expect_equal(centers$center_tract_id[centers$zone_id == "Z1"], "T2")
})

test_that(".find_zone_centers handles single-tract zones", {
  za <- surveyzones:::.validate_zone_assignments(data.frame(
    tract_id = "T1", zone_id = "Z1"
  ))
  centers <- surveyzones:::.find_zone_centers(za, tibble::tibble(
    origin_id = character(0), destination_id = character(0), distance = numeric(0)
  ))
  expect_equal(centers$center_tract_id, "T1")
  expect_equal(centers$n_tracts, 1L)
})

test_that(".find_zone_centers uses symmetrized (averaged) distances", {
  # A-B: only d(A,B)=1 given; A-C: d(A,C)=10 and d(C,A)=10;
  # B-C: d(B,C)=1 and d(C,B)=9 -> averaged 5
  # Sums: A = 1+10 = 11; B = 1+5 = 6; C = 10+5 = 15 -> medoid = B
  za <- surveyzones:::.validate_zone_assignments(data.frame(
    tract_id = c("A", "B", "C"), zone_id = c("Z1", "Z1", "Z1")
  ))
  dists <- tibble::tibble(
    origin_id      = c("A", "A", "C", "B", "C"),
    destination_id = c("B", "C", "A", "C", "B"),
    distance       = c(1, 10, 10, 1, 9)
  )
  centers <- surveyzones:::.find_zone_centers(za, dists)
  expect_equal(centers$center_tract_id, "B")
})

test_that(".find_zone_centers breaks ties deterministically (lowest tract_id)", {
  za <- surveyzones:::.validate_zone_assignments(data.frame(
    tract_id = c("T1", "T2", "T3"), zone_id = c("Z1", "Z1", "Z1")
  ))
  dists <- tibble::tibble(
    origin_id      = c("T1", "T1", "T2", "T2", "T3", "T3"),
    destination_id = c("T2", "T3", "T1", "T3", "T1", "T2"),
    distance       = c(2, 2, 2, 2, 2, 2)
  )
  centers <- surveyzones:::.find_zone_centers(za, dists)
  expect_equal(centers$center_tract_id, "T1")
})

test_that(".find_zone_centers aborts on incomplete within-zone distances", {
  za <- surveyzones:::.validate_zone_assignments(data.frame(
    tract_id = c("T1", "T2", "T3"), zone_id = c("Z1", "Z1", "Z1")
  ))
  # T1-T3 pair missing entirely
  dists <- tibble::tibble(
    origin_id      = c("T1", "T2"),
    destination_id = c("T2", "T1"),
    distance       = c(2, 2)
  )
  expect_error(
    surveyzones:::.find_zone_centers(za, dists),
    "incomplete"
  )
})

test_that(".find_zone_centers lists all incomplete zones in one error", {
  za <- surveyzones:::.validate_zone_assignments(data.frame(
    tract_id = c("T1", "T2", "T3", "T4"),
    zone_id  = c("Z1", "Z1", "Z2", "Z2")
  ))
  dists <- tibble::tibble(
    origin_id = character(0), destination_id = character(0), distance = numeric(0)
  )
  err <- tryCatch(
    surveyzones:::.find_zone_centers(za, dists),
    error = function(e) e
  )
  expect_true(grepl("Z1", conditionMessage(err)))
  expect_true(grepl("Z2", conditionMessage(err)))
})

test_that(".zone_dist_sparse melts the zone distance matrix without self-pairs", {
  assignments <- tibble::tibble(
    tract_id = paste0("T", 1:4),
    zone_id  = c("Z1", "Z1", "Z2", "Z2"),
    partition_id = "P1"
  )
  positions <- c(T1 = 0, T2 = 1, T3 = 2, T4 = 3)
  ids <- names(positions)
  dists <- tidyr::expand_grid(origin_id = ids, destination_id = ids) |>
    dplyr::filter(origin_id != destination_id) |>
    dplyr::mutate(distance = abs(positions[origin_id] - positions[destination_id]))

  out <- surveyzones:::.zone_dist_sparse(c("Z1", "Z2"), assignments, dists)

  expect_true(all(c("origin_id", "destination_id", "distance") %in% names(out)))
  expect_false(any(out$origin_id == out$destination_id))
  # min pairwise: T2-T3 = 1
  expect_equal(unique(out$distance), 1)
})

test_that(".zone_dist_sparse drops unreachable (Inf) zone pairs", {
  assignments <- tibble::tibble(
    tract_id = c("T1", "T2", "T3"),
    zone_id  = c("Z1", "Z2", "Z3"),
    partition_id = "P1"
  )
  # No distance data involving T3 at all
  dists <- tibble::tibble(
    origin_id = c("T1", "T2"), destination_id = c("T2", "T1"), distance = c(5, 5)
  )
  out <- surveyzones:::.zone_dist_sparse(c("Z1", "Z2", "Z3"), assignments, dists)
  expect_false(any(out$origin_id == "Z3" | out$destination_id == "Z3"))
})

test_that(".group_zones_pmedian binds close zones into one group", {
  skip_if_not_installed("ROI.plugin.glpk")
  za <- surveyzones:::.validate_zone_assignments(data.frame(
    tract_id = paste0("T", 1:4),
    zone_id  = paste0("Z", 1:4)
  ))
  positions <- c(T1 = 0, T2 = 1, T3 = 2, T4 = 3)
  ids <- names(positions)
  dists <- tidyr::expand_grid(origin_id = ids, destination_id = ids) |>
    dplyr::filter(origin_id != destination_id) |>
    dplyr::mutate(distance = abs(positions[origin_id] - positions[destination_id]))

  out <- surveyzones:::.group_zones_pmedian(
    zone_ids = paste0("Z", 1:4), zone_assignments = za, sparse_distances = dists,
    D_max = 10, solver = "glpk", max_time = 30, rel_tol = 0.01,
    verbose = FALSE, use_cache = FALSE
  )

  expect_equal(nrow(out$groups), 1L)
  expect_equal(out$groups$n_zones, 4L)
})

test_that(".group_zones_pmedian keeps every zone within D_max of its group center", {
  skip_if_not_installed("ROI.plugin.glpk")
  set.seed(1)
  n <- 8
  positions <- stats::runif(n, 0, 50)
  ids <- paste0("Z", seq_len(n))
  names(positions) <- ids
  za <- surveyzones:::.validate_zone_assignments(data.frame(
    tract_id = ids, zone_id = ids
  ))
  dists <- tidyr::expand_grid(origin_id = ids, destination_id = ids) |>
    dplyr::filter(origin_id != destination_id) |>
    dplyr::mutate(distance = abs(positions[origin_id] - positions[destination_id]))

  out <- surveyzones:::.group_zones_pmedian(
    zone_ids = ids, zone_assignments = za, sparse_distances = dists,
    D_max = 10, solver = "glpk", max_time = 30, rel_tol = 0.01,
    verbose = FALSE, use_cache = FALSE
  )

  dist_lookup <- dists |>
    dplyr::mutate(key = paste(origin_id, destination_id)) |>
    dplyr::select(key, distance) |>
    tibble::deframe()

  checked <- out$assignments |>
    dplyr::left_join(out$groups, by = "group_id") |>
    dplyr::filter(.data$zone_id != .data$center_zone_id) |>
    dplyr::mutate(d = dist_lookup[paste(.data$zone_id, .data$center_zone_id)])

  expect_true(all(checked$d <= 10))
})

test_that(".group_zones_pmedian isolates an out-of-range zone as a singleton", {
  skip_if_not_installed("ROI.plugin.glpk")
  ids <- c("Z1", "Z2", "Z3", "Z4")
  positions <- c(Z1 = 0, Z2 = 1, Z3 = 2, Z4 = 1000)
  za <- surveyzones:::.validate_zone_assignments(data.frame(
    tract_id = ids, zone_id = ids
  ))
  dists <- tidyr::expand_grid(origin_id = ids, destination_id = ids) |>
    dplyr::filter(origin_id != destination_id) |>
    dplyr::mutate(distance = abs(positions[origin_id] - positions[destination_id]))

  out <- surveyzones:::.group_zones_pmedian(
    zone_ids = ids, zone_assignments = za, sparse_distances = dists,
    D_max = 10, solver = "glpk", max_time = 30, rel_tol = 0.01,
    verbose = FALSE, use_cache = FALSE
  )

  z4_group <- out$assignments$group_id[out$assignments$zone_id == "Z4"]
  n_zones_in_z4_group <- out$groups$n_zones[out$groups$group_id == z4_group]
  expect_equal(n_zones_in_z4_group, 1L)
})

test_that(".group_zones_pmedian handles a single-zone partition trivially", {
  za <- surveyzones:::.validate_zone_assignments(data.frame(
    tract_id = "T1", zone_id = "Z1"
  ))
  out <- surveyzones:::.group_zones_pmedian(
    zone_ids = "Z1", zone_assignments = za,
    sparse_distances = tibble::tibble(
      origin_id = character(0), destination_id = character(0), distance = numeric(0)
    ),
    D_max = 10, solver = "glpk", max_time = 30, rel_tol = 0.01,
    verbose = FALSE, use_cache = FALSE
  )
  expect_equal(nrow(out$groups), 1L)
  expect_equal(out$groups$group_id, "Z1")
})

test_that(".group_zones_pmedian with solver='spopt' matches glpk's grouping structure", {
  skip_if_not_installed("ROI.plugin.glpk")
  skip_if_not_installed("spopt")
  za <- surveyzones:::.validate_zone_assignments(data.frame(
    tract_id = paste0("T", 1:4),
    zone_id  = paste0("Z", 1:4)
  ))
  positions <- c(T1 = 0, T2 = 1, T3 = 2, T4 = 3)
  ids <- names(positions)
  dists <- tidyr::expand_grid(origin_id = ids, destination_id = ids) |>
    dplyr::filter(origin_id != destination_id) |>
    dplyr::mutate(distance = abs(positions[origin_id] - positions[destination_id]))

  out_glpk <- surveyzones:::.group_zones_pmedian(
    zone_ids = paste0("Z", 1:4), zone_assignments = za, sparse_distances = dists,
    D_max = 10, solver = "glpk", max_time = 30, rel_tol = 0.01,
    verbose = FALSE, use_cache = FALSE
  )
  out_spopt <- surveyzones:::.group_zones_pmedian(
    zone_ids = paste0("Z", 1:4), zone_assignments = za, sparse_distances = dists,
    D_max = 10, solver = "spopt", max_time = 30, rel_tol = 0.01,
    verbose = FALSE, use_cache = FALSE
  )

  expect_equal(nrow(out_glpk$groups), nrow(out_spopt$groups))
  expect_equal(out_glpk$groups$n_zones, out_spopt$groups$n_zones)
})

test_that(".group_zones_pmedian use_cache=TRUE returns identical results on repeat calls", {
  skip_if_not_installed("ROI.plugin.glpk")
  za <- surveyzones:::.validate_zone_assignments(data.frame(
    tract_id = paste0("T", 1:4),
    zone_id  = paste0("Z", 1:4)
  ))
  positions <- c(T1 = 0, T2 = 1, T3 = 2, T4 = 3)
  ids <- names(positions)
  dists <- tidyr::expand_grid(origin_id = ids, destination_id = ids) |>
    dplyr::filter(origin_id != destination_id) |>
    dplyr::mutate(distance = abs(positions[origin_id] - positions[destination_id]))

  args <- list(
    zone_ids = paste0("Z", 1:4), zone_assignments = za, sparse_distances = dists,
    D_max = 10, solver = "glpk", max_time = 30, rel_tol = 0.01,
    verbose = FALSE, use_cache = TRUE
  )
  first <- do.call(surveyzones:::.group_zones_pmedian, args)
  second <- do.call(surveyzones:::.group_zones_pmedian, args)

  expect_equal(first$groups, second$groups)
  expect_equal(first$assignments, second$assignments)
})

test_that(".label_zone_groups produces 3-digit zzz and ggg permutations", {
  ids <- paste0("Z", 1:4)
  assignments <- tibble::tibble(
    tract_id = paste0("T", 1:4), zone_id = ids, partition_id = "P1"
  )
  positions <- c(T1 = 0, T2 = 1, T3 = 2, T4 = 3)
  dists <- tidyr::expand_grid(origin_id = names(positions), destination_id = names(positions)) |>
    dplyr::filter(origin_id != destination_id) |>
    dplyr::mutate(distance = abs(positions[origin_id] - positions[destination_id]))

  groups <- tibble::tibble(
    group_id = c("G1", "G2"),
    center_zone_id = c("Z1", "Z4"),
    n_zones = c(2L, 2L), diameter = c(1, 1)
  )

  out <- surveyzones:::.label_zone_groups(ids, groups, assignments, dists)

  expect_equal(sort(out$zone_labels$zzz), sprintf("%03d", 1:4))
  expect_equal(sort(out$group_labels$ggg), sprintf("%03d", 1:2))
})

test_that(".label_zone_groups keeps zzz independent of group membership", {
  # 4 zones on a line: Z1(0) Z2(1) Z3(2) Z4(3)
  # Groups: G1 = {Z1, Z3}, G2 = {Z2, Z4} -- interleaved on the line
  ids <- paste0("Z", 1:4)
  assignments <- tibble::tibble(tract_id = ids, zone_id = ids, partition_id = "P1")
  positions <- c(Z1 = 0, Z2 = 1, Z3 = 2, Z4 = 3)
  dists <- tidyr::expand_grid(origin_id = ids, destination_id = ids) |>
    dplyr::filter(origin_id != destination_id) |>
    dplyr::mutate(distance = abs(positions[origin_id] - positions[destination_id]))

  groups <- tibble::tibble(
    group_id = c("G1", "G2"), center_zone_id = c("Z1", "Z2"),
    n_zones = c(2L, 2L), diameter = c(2, 2)
  )

  out <- surveyzones:::.label_zone_groups(ids, groups, assignments, dists)

  zzz_num <- stats::setNames(as.integer(out$zone_labels$zzz), out$zone_labels$zone_id)

  # The zone TSP recovers the spatial line order (Z1..Z4 or its reverse)
  ordered_by_zzz <- names(sort(zzz_num))
  in_position_order <- all(ordered_by_zzz == paste0("Z", 1:4)) ||
    all(ordered_by_zzz == paste0("Z", 4:1))
  expect_true(in_position_order)

  # G1 = {Z1, Z3}: their zzz positions are 2 apart, not adjacent -- proving
  # the zone TSP was not constrained by group membership
  expect_equal(abs(zzz_num[["Z1"]] - zzz_num[["Z3"]]), 2L)
})

test_that(".label_zone_groups handles a single zone and single group", {
  ids <- "Z1"
  assignments <- tibble::tibble(tract_id = "T1", zone_id = "Z1", partition_id = "P1")
  dists <- tibble::tibble(
    origin_id = character(0), destination_id = character(0), distance = numeric(0)
  )
  groups <- tibble::tibble(
    group_id = "Z1", center_zone_id = "Z1", n_zones = 1L, diameter = 0
  )

  out <- surveyzones:::.label_zone_groups(ids, groups, assignments, dists)
  expect_equal(out$zone_labels$zzz, "001")
  expect_equal(out$group_labels$ggg, "001")
})

test_that(".anchor_partition_order reverses order when the last element is closer to start", {
  ids <- c("Z1", "Z2", "Z3")
  zone_centers <- c(Z1 = "T1", Z2 = "T2", Z3 = "T3")
  access_points <- sf::st_as_sf(
    data.frame(
      tract_id = c("T1", "T2", "T3"),
      lon = c(0, 1, 2), lat = c(0, 0, 0)
    ),
    coords = c("lon", "lat"), crs = 4326
  )
  starting_point <- sf::st_as_sf(
    data.frame(lon = 3, lat = 0), coords = c("lon", "lat"), crs = 4326
  )

  out <- surveyzones:::.anchor_partition_order(ids, zone_centers, starting_point, access_points)
  expect_equal(out, c("Z3", "Z2", "Z1"))
})

test_that(".anchor_partition_order keeps order when the first element is already closer", {
  ids <- c("Z1", "Z2", "Z3")
  zone_centers <- c(Z1 = "T1", Z2 = "T2", Z3 = "T3")
  access_points <- sf::st_as_sf(
    data.frame(
      tract_id = c("T1", "T2", "T3"),
      lon = c(0, 1, 2), lat = c(0, 0, 0)
    ),
    coords = c("lon", "lat"), crs = 4326
  )
  starting_point <- sf::st_as_sf(
    data.frame(lon = -1, lat = 0), coords = c("lon", "lat"), crs = 4326
  )

  out <- surveyzones:::.anchor_partition_order(ids, zone_centers, starting_point, access_points)
  expect_equal(out, c("Z1", "Z2", "Z3"))
})

test_that(".anchor_partition_order is a no-op for a single element or NULL starting_point", {
  expect_equal(
    surveyzones:::.anchor_partition_order("Z1", c(Z1 = "T1"), NULL, NULL),
    "Z1"
  )
  expect_equal(
    surveyzones:::.anchor_partition_order(c("Z1", "Z2"), c(Z1 = "T1", Z2 = "T2"), NULL, NULL),
    c("Z1", "Z2")
  )
})

test_that(".label_zone_groups anchors both zzz and ggg to the partition's starting point", {
  # 4 zones on a line: Z1(0) Z2(1) Z3(2) Z4(3); groups G1={Z1,Z2}, G2={Z3,Z4}.
  # Starting point is right at Z4 -> zzz should start at Z4 and ggg should
  # start with G2 (whose center is Z4).
  ids <- paste0("Z", 1:4)
  assignments <- tibble::tibble(tract_id = ids, zone_id = ids, partition_id = "P1")
  positions <- c(Z1 = 0, Z2 = 1, Z3 = 2, Z4 = 3)
  dists <- tidyr::expand_grid(origin_id = ids, destination_id = ids) |>
    dplyr::filter(origin_id != destination_id) |>
    dplyr::mutate(distance = abs(positions[origin_id] - positions[destination_id]))

  groups <- tibble::tibble(
    group_id = c("G1", "G2"), center_zone_id = c("Z1", "Z4"),
    n_zones = c(2L, 2L), diameter = c(1, 1)
  )

  zone_centers <- c(Z1 = "Z1", Z2 = "Z2", Z3 = "Z3", Z4 = "Z4")
  access_points <- sf::st_as_sf(
    data.frame(tract_id = ids, lon = unname(positions), lat = 0),
    coords = c("lon", "lat"), crs = 4326
  )
  starting_point <- sf::st_as_sf(
    data.frame(lon = 3, lat = 0), coords = c("lon", "lat"), crs = 4326
  )

  out <- surveyzones:::.label_zone_groups(
    ids, groups, assignments, dists,
    zone_centers = zone_centers, starting_point = starting_point,
    access_points = access_points
  )

  expect_equal(out$zone_labels$zone_id[out$zone_labels$zzz == "001"], "Z4")
  expect_equal(out$group_labels$group_id[out$group_labels$ggg == "001"], "G2")
})

test_that(".sequence_group_tracts produces a permutation of each zone's tracts", {
  za <- tibble::tibble(
    tract_id = paste0("T", 1:4),
    zone_id  = c("ZA", "ZA", "ZB", "ZB")
  )
  positions <- c(T1 = 0, T2 = 1, T3 = 10, T4 = 11)
  ids <- names(positions)
  dists <- tidyr::expand_grid(origin_id = ids, destination_id = ids) |>
    dplyr::filter(origin_id != destination_id) |>
    dplyr::mutate(distance = abs(positions[origin_id] - positions[destination_id]))

  out <- surveyzones:::.sequence_group_tracts(c("ZA", "ZB"), za, dists)

  for (zid in c("ZA", "ZB")) {
    assigned <- sort(za$tract_id[za$zone_id == zid])
    sequenced <- sort(out$tract_id[out$zone_id == zid])
    expect_equal(assigned, sequenced)
  }
})

test_that(".sequence_group_tracts orients the second zone toward the first zone's exit", {
  za <- tibble::tibble(
    tract_id = paste0("T", 1:6),
    zone_id  = c("ZA", "ZA", "ZA", "ZB", "ZB", "ZB")
  )
  positions <- c(T1 = 0, T2 = 1, T3 = 2, T4 = 10, T5 = 11, T6 = 12)
  ids <- names(positions)
  dists <- tidyr::expand_grid(origin_id = ids, destination_id = ids) |>
    dplyr::filter(origin_id != destination_id) |>
    dplyr::mutate(distance = abs(positions[origin_id] - positions[destination_id]))

  out <- surveyzones:::.sequence_group_tracts(c("ZA", "ZB"), za, dists)

  za_tracts <- out |> dplyr::filter(zone_id == "ZA") |> dplyr::arrange(visit_order)
  zb_tracts <- out |> dplyr::filter(zone_id == "ZB") |> dplyr::arrange(visit_order)

  exit_pos <- positions[za_tracts$tract_id[nrow(za_tracts)]]
  entry_pos <- positions[zb_tracts$tract_id[1]]
  tail_pos <- positions[zb_tracts$tract_id[nrow(zb_tracts)]]

  expect_true(abs(entry_pos - exit_pos) <= abs(tail_pos - exit_pos))
})

test_that(".sequence_group_tracts handles single-tract zones", {
  za <- tibble::tibble(tract_id = c("T1", "T2"), zone_id = c("Z1", "Z2"))
  dists <- tibble::tibble(
    origin_id = character(0), destination_id = character(0), distance = numeric(0)
  )
  out <- surveyzones:::.sequence_group_tracts(c("Z1", "Z2"), za, dists)
  expect_equal(nrow(out), 2L)
  expect_equal(out$visit_order, c(1L, 1L))
})

test_that("new_surveyzones_zone_groups and print work", {
  obj <- surveyzones:::new_surveyzones_zone_groups(
    zones = tibble::tibble(
      zone_id = "Z1", partition_id = "P1", center_tract_id = "T1",
      n_tracts = 1L, group_id = "Z1", ggg = "001", zzz = "001", label = "001.001"
    ),
    groups = tibble::tibble(
      group_id = "Z1", partition_id = "P1", center_zone_id = "Z1",
      n_zones = 1L, diameter = 0, ggg = "001"
    ),
    sequence = tibble::tibble(zone_id = "Z1", tract_id = "T1", visit_order = 1L),
    parameters = list(D_max = 10, solver = "glpk", timestamp = Sys.time()),
    diagnostics = list(solver_status = "optimal")
  )
  expect_s3_class(obj, "surveyzones_zone_groups")
  expect_no_error(print(obj))
})

test_that("surveyzones_group_zones returns a surveyzones_zone_groups object", {
  skip_if_not_installed("ROI.plugin.glpk")
  za <- data.frame(
    tract_id = paste0("T", 1:4),
    zone_id  = c("Z1", "Z1", "Z2", "Z2")
  )
  positions <- c(T1 = 0, T2 = 1, T3 = 2, T4 = 3)
  ids <- names(positions)
  dists <- tidyr::expand_grid(origin_id = ids, destination_id = ids) |>
    dplyr::filter(origin_id != destination_id) |>
    dplyr::mutate(distance = abs(positions[origin_id] - positions[destination_id]))

  result <- surveyzones_group_zones(
    zone_assignments = za,
    sparse_distances = dists,
    D_max = 10,
    solver = "glpk",
    use_cache = FALSE,
    complete_distances = FALSE
  )

  expect_s3_class(result, "surveyzones_zone_groups")
  expect_equal(nrow(result$zones), 2L)
  expect_equal(sort(result$sequence$tract_id), sort(za$tract_id))
})

test_that("surveyzones_group_zones requires access_points when complete_distances = TRUE", {
  za <- data.frame(tract_id = "T1", zone_id = "Z1")
  dists <- tibble::tibble(
    origin_id = character(0), destination_id = character(0), distance = numeric(0)
  )
  expect_error(
    surveyzones_group_zones(za, dists, D_max = 10),
    "access_points"
  )
})

test_that("surveyzones_group_zones resolves centers after haversine completion", {
  skip_if_not_installed("ROI.plugin.glpk")
  za <- data.frame(
    tract_id = c("T1", "T2", "T3"),
    zone_id  = c("Z1", "Z1", "Z1")
  )
  pts <- sf::st_as_sf(
    data.frame(
      tract_id = c("T1", "T2", "T3"),
      lon = c(-38.50, -38.501, -38.502),
      lat = c(-13.00, -13.001, -13.002)
    ),
    coords = c("lon", "lat"), crs = 4326
  )
  # Only one pair supplied; the rest must come from haversine completion
  dists <- tibble::tibble(
    origin_id = c("T1", "T2"), destination_id = c("T2", "T1"), distance = c(1, 1)
  )

  result <- surveyzones_group_zones(
    zone_assignments = za,
    sparse_distances = dists,
    D_max = 100,
    solver = "glpk",
    use_cache = FALSE,
    access_points = pts,
    complete_distances = TRUE
  )

  expect_equal(nrow(result$zones), 1L)
  expect_true(result$zones$center_tract_id %in% za$tract_id)
})

test_that("surveyzones_group_zones anchors the first zone/group per partition to starting_points", {
  skip_if_not_installed("ROI.plugin.glpk")
  # Two groups on a line: {Z1,Z2} near position 0-1, {Z3,Z4} near position
  # 10-11. D_max=5 keeps them as two separate groups (Z2-Z3 distance = 9 > 5).
  # Starting point at position 11 (right at Z4) -> Z4 must be the first zone
  # (zzz == "001"), and the group CONTAINING Z4 must be first (ggg == "001").
  # We deliberately do not assert which zone becomes center_zone_id within
  # that group -- anchoring reorders the group/zone sequences, it does not
  # change which zone the p-median solver picked as a group's center.
  za <- tibble::tibble(
    zone_id = paste0("Z", 1:4),
    tract_id = paste0("Z", 1:4),
    partition_id = "P1"
  )
  positions <- c(Z1 = 0, Z2 = 1, Z3 = 10, Z4 = 11)
  pts <- sf::st_as_sf(
    data.frame(tract_id = names(positions), lon = unname(positions), lat = 0),
    coords = c("lon", "lat"), crs = 4326
  )
  dists <- tidyr::expand_grid(origin_id = names(positions), destination_id = names(positions)) |>
    dplyr::filter(origin_id != destination_id) |>
    dplyr::mutate(distance = abs(positions[origin_id] - positions[destination_id]))

  starting_points <- sf::st_as_sf(
    data.frame(partition_id = "P1", lon = 11, lat = 0),
    coords = c("lon", "lat"), crs = 4326
  )

  result <- surveyzones_group_zones(
    zone_assignments = za,
    sparse_distances = dists,
    D_max = 5,
    solver = "glpk",
    use_cache = FALSE,
    access_points = pts,
    complete_distances = FALSE,
    starting_points = starting_points
  )

  expect_equal(result$zones$zone_id[result$zones$zzz == "001"], "Z4")

  z4_group <- result$zones$group_id[result$zones$zone_id == "Z4"]
  expect_equal(unique(result$groups$ggg[result$groups$group_id == z4_group]), "001")
})

test_that("surveyzones_group_zones end-to-end: shapes, D_max bound, labels, tract coverage", {
  skip_if_not_installed("ROI.plugin.glpk")

  # 7 zones across 2 partitions, positions on a line within each partition.
  # Partition A: zones at 0,1,2,10,11 (5 zones); D_max=5 groups {0,1,2} and {10,11}.
  # Partition B: zones at 0,50 (2 zones); D_max=5 -> 2 singleton groups.
  zone_positions <- c(
    A_Z1 = 0, A_Z2 = 1, A_Z3 = 2, A_Z4 = 10, A_Z5 = 11,
    B_Z1 = 0, B_Z2 = 50
  )
  zone_partition <- c(rep("A", 5), rep("B", 2))
  names(zone_partition) <- names(zone_positions)

  # 2 tracts per zone, both located at the zone's position
  za <- tibble::tibble(
    zone_id = rep(names(zone_positions), each = 2),
    tract_id = paste0(rep(names(zone_positions), each = 2), "_t", 1:2),
    partition_id = rep(zone_partition, each = 2)
  )
  tract_pos <- stats::setNames(rep(zone_positions, each = 2), za$tract_id)

  ids <- za$tract_id
  dists <- tidyr::expand_grid(origin_id = ids, destination_id = ids) |>
    dplyr::filter(origin_id != destination_id) |>
    dplyr::mutate(distance = abs(tract_pos[origin_id] - tract_pos[destination_id]))

  result <- surveyzones_group_zones(
    zone_assignments = za,
    sparse_distances = dists,
    D_max = 5,
    solver = "glpk",
    use_cache = FALSE,
    complete_distances = FALSE
  )

  expect_s3_class(result, "surveyzones_zone_groups")
  expect_true(all(c("zone_id", "partition_id", "center_tract_id", "n_tracts",
                     "group_id", "ggg", "zzz", "label") %in% names(result$zones)))
  expect_true(all(c("group_id", "partition_id", "center_zone_id", "n_zones",
                     "diameter", "ggg") %in% names(result$groups)))
  expect_true(all(c("zone_id", "tract_id", "visit_order") %in% names(result$sequence)))

  expect_equal(nrow(result$zones), 7L)
  expect_equal(sort(result$sequence$tract_id), sort(za$tract_id))

  # Label format
  expect_true(all(grepl("^\\d{3}\\.\\d{3}$", result$zones$label)))

  # ggg and zzz restart at 001 in each partition
  for (pid in c("A", "B")) {
    part_zones <- result$zones[result$zones$partition_id == pid, ]
    expect_true("001" %in% part_zones$zzz)
    expect_true("001" %in% part_zones$ggg)
  }

  # Partition A groups into 2 clusters (close triple + close pair)
  a_groups <- unique(result$zones$group_id[result$zones$partition_id == "A"])
  expect_equal(length(a_groups), 2L)

  # Partition B: two zones far apart -> 2 singleton groups
  b_groups <- unique(result$zones$group_id[result$zones$partition_id == "B"])
  expect_equal(length(b_groups), 2L)

  # No group spans two partitions: every zone in a given group_id has the
  # same partition_id
  partitions_per_group <- result$zones |>
    dplyr::summarise(n_partitions = dplyr::n_distinct(.data$partition_id), .by = "group_id")
  expect_true(all(partitions_per_group$n_partitions == 1L))
})

test_that("surveyzones_group_zones is deterministic across repeat calls", {
  skip_if_not_installed("ROI.plugin.glpk")
  # NOTE: without a starting_points anchor, an open-path TSP has an inherent
  # forward/reverse tie (both directions have identical total length), so
  # .seriate_1d()'s underlying heuristic can legitimately return either
  # direction on different calls -- this is not a bug (the same
  # non-determinism was independently confirmed for surveyzones_sequence_tracts()
  # in Task 6's review). Anchoring with starting_points breaks that tie
  # deterministically (verified empirically: 3 repeat calls with the same
  # starting_points produced byte-identical zzz/ggg/group_id every time),
  # so this test anchors both partitions to make the full pipeline's output
  # deterministic, which is the property actually guaranteed by the design.
  za <- tibble::tibble(
    zone_id = paste0("Z", 1:6),
    tract_id = paste0("Z", 1:6),
    partition_id = c(rep("A", 4), rep("B", 2))
  )
  positions <- c(Z1 = 0, Z2 = 1, Z3 = 2, Z4 = 3, Z5 = 0, Z6 = 20)
  ids <- names(positions)
  pts <- sf::st_as_sf(
    data.frame(tract_id = ids, lon = unname(positions), lat = 0),
    coords = c("lon", "lat"), crs = 4326
  )
  dists <- tidyr::expand_grid(origin_id = ids, destination_id = ids) |>
    dplyr::filter(origin_id != destination_id) |>
    dplyr::mutate(distance = abs(positions[origin_id] - positions[destination_id]))

  starting_points <- sf::st_as_sf(
    data.frame(partition_id = c("A", "B"), lon = c(1.5, 10), lat = 0),
    coords = c("lon", "lat"), crs = 4326
  )

  args <- list(
    zone_assignments = za, sparse_distances = dists, D_max = 5,
    solver = "glpk", use_cache = FALSE, access_points = pts,
    complete_distances = FALSE, starting_points = starting_points
  )
  set.seed(42)
  first <- do.call(surveyzones_group_zones, args)
  set.seed(42)
  second <- do.call(surveyzones_group_zones, args)

  expect_equal(first$zones, second$zones)
  expect_equal(first$groups, second$groups)
  expect_equal(first$sequence, second$sequence)
})

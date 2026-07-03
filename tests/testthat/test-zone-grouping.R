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

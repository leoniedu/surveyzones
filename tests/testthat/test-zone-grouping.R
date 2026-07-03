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

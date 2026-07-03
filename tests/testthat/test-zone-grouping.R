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

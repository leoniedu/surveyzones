test_that("partitioning splits by partition_id", {
  tracts <- data.frame(
    tract_id = c("A", "B", "C", "D"),
    partition_id = c("P1", "P1", "P2", "P2"),
    expected_service_time = c(1, 1, 1, 1)
  )

  parts <- surveyzones_partition(tracts, enforce_partition = TRUE)

  expect_type(parts, "list")
  expect_equal(length(parts), 2)
  expect_true(all(c("P1", "P2") %in% names(parts)))
  expect_equal(nrow(parts[["P1"]]), 2)
  expect_equal(nrow(parts[["P2"]]), 2)
})

test_that("partitioning returns single group when no partition_id", {
  tracts <- data.frame(
    tract_id = c("A", "B", "C"),
    expected_service_time = c(1, 1, 1)
  )

  parts <- surveyzones_partition(tracts, enforce_partition = TRUE)

  expect_equal(length(parts), 1)
  expect_equal(names(parts), "all")
  expect_equal(nrow(parts[["all"]]), 3)
})

test_that("enforce_partition = FALSE ignores partition_id", {
  tracts <- data.frame(
    tract_id = c("A", "B", "C"),
    partition_id = c("P1", "P1", "P2"),
    expected_service_time = c(1, 1, 1)
  )

  parts <- surveyzones_partition(tracts, enforce_partition = FALSE)

  expect_equal(length(parts), 1)
  expect_equal(names(parts), "all")
})

test_that("partitioning rejects missing tract_id", {
  tracts <- data.frame(id = c("A", "B"))
  expect_error(surveyzones_partition(tracts), "tract_id")
})

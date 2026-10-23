test_that("mandatory tracts land in roteiros, at most k per team, band respected", {
  skip_if_not_installed("ROI.plugin.highs")
  x <- toy_instance()
  res <- surveyzones_assign_teams(x$units, x$points, x$teams, x$eligible, x$durations,
                                  mandatory = c("A", "B", "E"), k = 2L, tol = 0.25, solver = "highs")
  expect_s3_class(res, "surveyzones_assignment")
  expect_setequal(res$assignments$tract_id, x$units$tract_id)
  expect_true(all(res$assignments$mode[res$assignments$tract_id %in% c("A", "B", "E")] == "roteiro"))
  expect_true(all(res$loads$n_roteiros <= 2L))
  expect_true(all(res$loads$workload >= res$loads$lower & res$loads$workload <= res$loads$upper))
})

test_that("ineligible pairs are never used and fixed pairs are kept", {
  skip_if_not_installed("ROI.plugin.highs")
  x <- toy_instance()
  elig <- dplyr::filter(x$eligible, !(tract_id == "A" & team_id == "t1"))
  res <- surveyzones_assign_teams(x$units, x$points, x$teams, elig, x$durations,
                                  fixed = tibble::tibble(tract_id = "F", team_id = "t1"),
                                  tol = 0.5, solver = "highs")
  a <- res$assignments
  expect_equal(a$team_id[a$tract_id == "A"], "t2")
  expect_equal(a$mode[a$tract_id == "F"], "fixed")
  expect_equal(a$team_id[a$tract_id == "F"], "t1")
})

test_that("optimal on a tiny instance (brute force)", {
  skip_if_not_installed("ROI.plugin.highs")
  x <- toy_instance(n = 4)
  res <- surveyzones_assign_teams(x$units, x$points, x$teams, x$eligible, x$durations,
                                  mandatory = x$units$tract_id, k = 2L, tol = 0.5, solver = "highs")
  expect_equal(res$diagnostics$objective_value, brute_force_objective(x, k = 2L, tol = 0.5), tolerance = 1e-6)
})

test_that("asymmetric durations are used origin -> centre", {
  skip_if_not_installed("ROI.plugin.highs")
  x <- toy_instance()
  res <- surveyzones_assign_teams(x$units, x$points, x$teams, x$eligible, x$durations, tol = 0.5, solver = "highs")
  r <- dplyr::left_join(res$roteiros, res$origins, by = dplyr::join_by(team_id))
  d <- dplyr::left_join(r, x$durations, by = dplyr::join_by(origin_tract_id == origin_id, center_tract_id == destination_id))
  expect_equal(r$cost, dplyr::coalesce(d$distance, 0))
})

test_that("same inputs give the same assignment", {
  skip_if_not_installed("ROI.plugin.highs")
  x <- toy_instance()
  r1 <- surveyzones_assign_teams(x$units, x$points, x$teams, x$eligible, x$durations, tol = 0.5, solver = "highs")
  r2 <- surveyzones_assign_teams(x$units, x$points, x$teams, x$eligible, x$durations, tol = 0.5, solver = "highs")
  expect_identical(r1$assignments, r2$assignments)
})

test_that("units missing from durations abort with their ids", {
  x <- toy_instance()
  expect_error(
    surveyzones_assign_teams(x$units, x$points, x$teams, x$eligible,
                             dplyr::filter(x$durations, origin_id != "C", destination_id != "C"), solver = "highs"),
    "C"
  )
})

test_that("infeasible hard solve returns an empty assignment", {
  skip_if_not_installed("ROI.plugin.highs")
  x <- toy_instance()
  res <- surveyzones_assign_teams(x$units, x$points, x$teams, x$eligible, x$durations,
                                  mandatory = x$units$tract_id, k = 1L, tol = 0.01, solver = "highs")
  expect_equal(res$diagnostics$solver_status, "infeasible")
  expect_equal(nrow(res$assignments), 0L)
})

test_that("invalid inputs abort", {
  x <- toy_instance()
  expect_error(surveyzones_assign_teams(x$units, x$points, x$teams, x$eligible, x$durations, k = 0L, solver = "highs"), "k")
  expect_error(surveyzones_assign_teams(x$units, x$points, x$teams, x$eligible, x$durations, tol = -1, solver = "highs"), "tol")
  expect_error(surveyzones_assign_teams(x$units, x$points, x$teams, x$eligible, x$durations, mandatory = "Z", solver = "highs"), "Z")
})

test_that("objective equals roteiro plus loose costs and loose tracts have anchors", {
  skip_if_not_installed("ROI.plugin.highs")
  x <- toy_instance()
  res <- surveyzones_assign_teams(x$units, x$points, x$teams, x$eligible, x$durations,
                                  mandatory = "A", tol = 0.5, max_tracts = 1L, solver = "highs")
  expect_equal(
    res$diagnostics$objective_value,
    sum(res$roteiros$cost) + sum(res$loads$loose_cost),
    tolerance = 1e-6
  )
  loose <- dplyr::filter(res$assignments, mode == "loose")
  expect_gt(nrow(loose), 0L)
  expect_true(all(loose$team_id %in% res$roteiros$team_id))
})

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
  res <- surveyzones_assign_teams(x$units, x$points, x$teams, x$eligible, x$durations,
                                  mandatory = x$units$tract_id, k = 3L, max_tracts = 1L, tol = 1, solver = "highs")
  r <- dplyr::left_join(res$roteiros, res$origins, by = dplyr::join_by(team_id))
  d <- dplyr::left_join(r, x$durations, by = dplyr::join_by(origin_tract_id == origin_id, center_tract_id == destination_id))
  expect_equal(r$cost, dplyr::coalesce(d$distance, 0))
  expect_true(any(r$cost > 0))
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

test_that("infeasible hard solve returns an elastic near-plan", {
  skip_if_not_installed("ROI.plugin.highs")
  x <- toy_instance()
  res <- surveyzones_assign_teams(x$units, x$points, x$teams, x$eligible, x$durations,
                                  mandatory = x$units$tract_id, k = 1L, tol = 0.01, solver = "highs")
  expect_equal(res$diagnostics$solver_status, "infeasible")
  expect_false(is.null(res$diagnostics$elastic))
  expect_gt(nrow(res$diagnostics$elastic), 0L)
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
  for (i in seq_len(nrow(loose))) {
    anchors <- res$assignments$tract_id[res$assignments$team_id == loose$team_id[i] & res$assignments$mode != "loose"]
    dur <- dplyr::filter(x$durations, destination_id == loose$tract_id[i], origin_id %in% anchors)
    expect_true(nrow(dur) > 0L)
    expect_equal(loose$cost[i], min(dur$distance), tolerance = 1e-6)
  }
})

test_that("floor shortfall is reported per team", {
  skip_if_not_installed("ROI.plugin.highs")
  x <- toy_instance()
  elig <- dplyr::filter(x$eligible, team_id == "t1" | tract_id == "A")
  res <- surveyzones_assign_teams(x$units, x$points, x$teams, elig, x$durations, tol = 0.1, solver = "highs")
  expect_equal(res$diagnostics$solver_status, "infeasible")
  e <- res$diagnostics$elastic
  expect_true(any(e$rule == "band_lower" & e$id == "t2"))
  expect_gt(sum(e$amount[e$rule == "band_lower" & e$id == "t2"]), 0)
})

test_that("feasible coverage but infeasible band is diagnosed, not crashed", {
  skip_if_not_installed("ROI.plugin.highs")
  x <- toy_instance()
  res <- surveyzones_assign_teams(x$units, x$points, x$teams, x$eligible, x$durations,
                                  mandatory = x$units$tract_id, k = 6L, tol = 0, solver = "highs")
  expect_equal(res$diagnostics$solver_status, "infeasible")
  expect_true(all(res$diagnostics$elastic$rule %in% c("band_lower", "band_upper")))
  expect_gt(nrow(res$assignments), 0L)
})

test_that("min_roteiros is the set-partition minimum", {
  skip_if_not_installed("ROI.plugin.highs")
  x <- toy_instance()
  res <- surveyzones_assign_teams(x$units, x$points, x$teams, x$eligible, x$durations,
                                  mandatory = c("A", "B", "C"), tol = 0.5, solver = "highs")
  expect_equal(res$diagnostics$min_roteiros, 1L)
})

test_that(".solve_op classifies infeasible and optimal via the solver status symbol", {
  skip_if_not_installed("ROI.plugin.highs")
  mk <- function(rhs) {
    ROI::OP(
      ROI::L_objective(c(1, 1)),
      ROI::L_constraint(matrix(c(1, 1), 1), "<=", rhs),
      types = c("B", "B"), maximum = TRUE
    )
  }
  bad <- .solve_op(mk(-1), "highs", max_time = 10, rel_tol = 0, verbose = FALSE)
  expect_equal(bad$status, "infeasible")
  expect_null(bad$solution)
  ok <- .solve_op(mk(1), "highs", max_time = 10, rel_tol = 0, verbose = FALSE)
  expect_equal(ok$status, "optimal")
  expect_equal(ok$objective_value, 1)
})

test_that("print shows the gap for feasible-but-not-optimal solves", {
  skip_if_not_installed("ROI.plugin.highs")
  x <- toy_instance()
  res <- surveyzones_assign_teams(x$units, x$points, x$teams, x$eligible, x$durations,
                                  k = 2L, tol = 0.5, solver = "highs")
  res$diagnostics$solver_status <- "feasible"
  res$diagnostics$optimality_gap <- 0.034
  expect_match(paste(cli::cli_fmt(print(res)), collapse = "\n"), "3.4%")
})

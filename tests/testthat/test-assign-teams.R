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

mock_roi_solution <- function(symbol = NULL, code = 1L, solution = c(1, 0), objval = 3,
                              message = NULL) {
  list(
    solution = solution, objval = objval,
    status = list(code = code, msg = list(symbol = symbol, message = "")),
    message = message
  )
}

test_that("ROI status classification handles HiGHS and glpk symbols", {
  cls <- function(...) .classify_roi_solution(mock_roi_solution(...))$status
  expect_equal(cls("OPTIMAL", code = 0L), "optimal")
  expect_equal(cls("GLP_OPT", code = 0L), "optimal")
  expect_equal(cls(NULL, code = 0L), "optimal")
  expect_equal(cls("TIME_LIMIT"), "feasible")
  expect_equal(cls("GLP_FEAS"), "feasible")
  for (s in c("INFEASIBLE", "UNBOUNDEDORINFEASIBLE", "GLP_NOFEAS", "GLP_INFEAS")) {
    expect_equal(cls(s, solution = NULL, objval = NA_real_), "infeasible")
  }
  expect_equal(cls("GLP_UNDEF", solution = NULL, objval = NA_real_), "error")
  expect_equal(cls("GLP_UNBND", solution = NULL, objval = NA_real_), "error")
  expect_equal(cls("GLP_OPT", code = 0L, solution = NULL, objval = NA_real_), "error")
  expect_equal(
    cls("WEIRD", solution = NULL, objval = NA_real_, message = list(status = "problem is infeasible")),
    "infeasible"
  )
})

test_that("k = 'auto' follows the minimum-roteiros formula", {
  skip_if_not_installed("ROI.plugin.highs")
  x <- toy_instance()
  res <- surveyzones_assign_teams(x$units, x$points, x$teams, x$eligible, x$durations,
                                  mandatory = x$units$tract_id, k = "auto", tol = 1, solver = "highs")
  k_min <- ceiling(res$diagnostics$min_roteiros / nrow(x$teams))
  expect_true(res$parameters$k_auto)
  expect_equal(res$parameters$k_min, k_min)
  expect_equal(res$parameters$k, k_min + 2L)
  expect_message(print(res), "auto; minimum")
})

test_that("integer k is unchanged and non-whole k aborts", {
  skip_if_not_installed("ROI.plugin.highs")
  x <- toy_instance()
  res <- surveyzones_assign_teams(x$units, x$points, x$teams, x$eligible, x$durations,
                                  k = 2L, tol = 0.5, solver = "highs")
  expect_equal(res$parameters$k, 2L)
  expect_false(res$parameters$k_auto)
  expect_error(
    surveyzones_assign_teams(x$units, x$points, x$teams, x$eligible, x$durations, k = 2.5, solver = "highs"),
    "whole number"
  )
})

local_assign_cache <- function(env = parent.frame()) {
  tmp <- withr::local_tempdir(.local_envir = env)
  other <- withr::local_tempdir(.local_envir = env)
  # .onLoad assigns the memoised function with `<<-`, which lands in globalenv
  old_mem <- get(".assign_teams_core_mem", envir = globalenv())
  assign(".assign_teams_core_mem", memoise::memoise(
    surveyzones:::.assign_teams_core,
    omit_args = "verbose",
    cache = cachem::cache_disk(dir = tmp)
  ), envir = globalenv())
  pkg <- surveyzones:::.pkg
  old <- as.list(pkg)[c("cache_dir", "assign_cache_dir")]
  pkg$cache_dir <- other
  pkg$assign_cache_dir <- tmp
  withr::defer({
    assign(".assign_teams_core_mem", old_mem, envir = globalenv())
    pkg$cache_dir <- old$cache_dir
    pkg$assign_cache_dir <- old$assign_cache_dir
  }, envir = env)
  tmp
}

test_that("use_cache: repeat call hits, changed input misses, use_cache = FALSE bypasses", {
  skip_if_not_installed("ROI.plugin.highs")
  local_assign_cache()
  x <- toy_instance()
  run <- function(...) {
    surveyzones_assign_teams(x$units, x$points, x$teams, x$eligible, x$durations,
                             k = 2L, solver = "highs", ...)
  }
  expect_message(r1 <- run(tol = 0.5), "Computing and caching")
  expect_message(r2 <- run(tol = 0.5), "Using cached")
  expect_identical(r1, r2)
  expect_message(r3 <- run(tol = 0.5, verbose = TRUE), "Using cached")
  expect_message(run(tol = 0.4), "Computing and caching")
  expect_no_message(r4 <- run(tol = 0.5, use_cache = FALSE), message = "cached")
  expect_equal(r4$assignments, r1$assignments)
  expect_equal(r4$diagnostics$objective_value, r1$diagnostics$objective_value)
})

test_that("surveyzones_clear_cache(force = TRUE) clears the assign_teams cache", {
  skip_if_not_installed("ROI.plugin.highs")
  tmp <- local_assign_cache()
  x <- toy_instance()
  run <- function() {
    surveyzones_assign_teams(x$units, x$points, x$teams, x$eligible, x$durations,
                             k = 2L, tol = 0.5, solver = "highs")
  }
  suppressMessages(run())
  expect_gt(length(list.files(tmp)), 0L)
  suppressMessages(surveyzones_clear_cache(force = TRUE))
  expect_length(list.files(tmp, recursive = TRUE), 0L)
  expect_message(run(), "Computing and caching")
})

.origin_cost <- function(x, res, origins) {
  vapply(seq_len(nrow(res$roteiros)), \(i) {
    og <- origins$tract_id[origins$team_id == res$roteiros$team_id[i]]
    ctr <- res$roteiros$center_tract_id[i]
    min(vapply(og, \(o) {
      if (o == ctr) return(0)
      x$durations$distance[x$durations$origin_id == o & x$durations$destination_id == ctr]
    }, numeric(1)))
  }, numeric(1))
}

test_that("given origins fix roteiro costs for all teams", {
  skip_if_not_installed("ROI.plugin.highs")
  x <- toy_instance()
  og <- tibble::tibble(team_id = c("t1", "t2"), tract_id = c("A", "F"))
  res <- surveyzones_assign_teams(x$units, x$points, x$teams, x$eligible, x$durations,
                                  mandatory = x$units$tract_id, origins = og, k = 3L, max_tracts = 1L,
                                  tol = 1, solver = "highs", use_cache = FALSE)
  expect_equal(res$roteiros$cost, .origin_cost(x, res, og))
  expect_true(any(res$roteiros$cost > 0))
  expect_equal(res$origins$source, c("given", "given"))
  expect_setequal(res$origins$origin_tract_id, c("A", "F"))
  expect_identical(res$parameters$origins, og)
})

test_that("mixed given and model origins both behave", {
  skip_if_not_installed("ROI.plugin.highs")
  x <- toy_instance()
  og <- tibble::tibble(team_id = "t1", tract_id = "A")
  res <- surveyzones_assign_teams(x$units, x$points, x$teams, x$eligible, x$durations,
                                  mandatory = x$units$tract_id, origins = og, k = 3L, max_tracts = 1L,
                                  tol = 1, solver = "highs", use_cache = FALSE)
  expect_equal(res$origins$source[res$origins$team_id == "t1"], "given")
  expect_equal(res$origins$source[res$origins$team_id == "t2"], "model")
  expect_equal(nrow(res$origins), 2L)
  r1 <- res$roteiros[res$roteiros$team_id == "t1", ]
  expect_equal(r1$cost, .origin_cost(x, list(roteiros = r1), og))
  r2 <- res$roteiros[res$roteiros$team_id == "t2", ]
  o2 <- res$origins$origin_tract_id[res$origins$team_id == "t2"]
  expect_equal(r2$cost, .origin_cost(x, list(roteiros = r2), tibble::tibble(team_id = "t2", tract_id = o2)))
  expect_equal(res$diagnostics$objective_value, sum(res$roteiros$cost) + sum(res$loads$loose_cost), tolerance = 1e-6)
})

test_that("two origins per team use the closer one", {
  skip_if_not_installed("ROI.plugin.highs")
  x <- toy_instance()
  og <- tibble::tibble(team_id = c("t1", "t1", "t2", "t2"), tract_id = c("A", "F", "B", "E"))
  res <- surveyzones_assign_teams(x$units, x$points, x$teams, x$eligible, x$durations,
                                  mandatory = x$units$tract_id, origins = og, k = 3L, max_tracts = 1L,
                                  tol = 1, solver = "highs", use_cache = FALSE)
  expect_equal(nrow(res$origins), 4L)
  expect_equal(res$roteiros$cost, .origin_cost(x, res, og))
})

test_that("origins outside units use their own durations", {
  skip_if_not_installed("ROI.plugin.highs")
  x <- toy_instance()
  casa <- tibble::tibble(
    origin_id = "casa1", destination_id = x$units$tract_id,
    distance = c(9, 9, 9, 9, 9, 0.25)
  )
  x$durations <- dplyr::bind_rows(x$durations, casa)
  og <- tibble::tibble(team_id = "t1", tract_id = "casa1")
  res <- surveyzones_assign_teams(x$units, x$points, x$teams, x$eligible, x$durations,
                                  mandatory = x$units$tract_id, origins = og, k = 3L, max_tracts = 1L,
                                  tol = 1, solver = "highs", use_cache = FALSE)
  r1 <- res$roteiros[res$roteiros$team_id == "t1", ]
  esperado <- casa$distance[match(r1$center_tract_id, casa$destination_id)]
  expect_equal(r1$cost, esperado)
  expect_equal(res$origins$origin_tract_id[res$origins$team_id == "t1"], "casa1")
  expect_equal(res$diagnostics$objective_value,
               sum(res$roteiros$cost) + sum(res$loads$loose_cost), tolerance = 1e-6)
})

test_that("optimal with given origins (brute force)", {
  skip_if_not_installed("ROI.plugin.highs")
  x <- toy_instance(n = 4)
  og <- tibble::tibble(team_id = c("t1", "t1", "t2"), tract_id = c("A", "D", "B"))
  res <- surveyzones_assign_teams(x$units, x$points, x$teams, x$eligible, x$durations,
                                  mandatory = x$units$tract_id, origins = og, k = 2L, tol = 0.5,
                                  solver = "highs", use_cache = FALSE)
  expect_equal(res$diagnostics$objective_value,
               brute_force_objective(x, k = 2L, tol = 0.5, origins = og), tolerance = 1e-6)
})

test_that("invalid origins abort naming offenders", {
  x <- toy_instance()
  run <- function(og) {
    surveyzones_assign_teams(x$units, x$points, x$teams, x$eligible, x$durations,
                             origins = og, solver = "highs", use_cache = FALSE)
  }
  expect_error(run(tibble::tibble(team_id = "tX", tract_id = "A")), "tX")
  expect_error(run(tibble::tibble(team_id = "t1", tract_id = "ZZ")), "ZZ")
  expect_error(run(tibble::tibble(team_id = "t1")), "origins")
})

test_that("origins work with the elastic fallback and min_roteiros", {
  skip_if_not_installed("ROI.plugin.highs")
  x <- toy_instance()
  og <- tibble::tibble(team_id = "t1", tract_id = "A")
  res <- surveyzones_assign_teams(x$units, x$points, x$teams, x$eligible, x$durations,
                                  mandatory = x$units$tract_id, origins = og, k = 1L, tol = 0.01,
                                  solver = "highs", use_cache = FALSE)
  expect_equal(res$diagnostics$solver_status, "infeasible")
  expect_gt(nrow(res$diagnostics$elastic), 0L)
  expect_true("given" %in% res$origins$source)
  expect_true(is.integer(res$diagnostics$min_roteiros))
  auto <- surveyzones_assign_teams(x$units, x$points, x$teams, x$eligible, x$durations,
                                   mandatory = x$units$tract_id, origins = og, k = "auto", tol = 1,
                                   solver = "highs", use_cache = FALSE)
  expect_true(auto$parameters$k_auto)
})

test_that("use_cache: changing origins misses", {
  skip_if_not_installed("ROI.plugin.highs")
  local_assign_cache()
  x <- toy_instance()
  run <- function(og) {
    surveyzones_assign_teams(x$units, x$points, x$teams, x$eligible, x$durations,
                             origins = og, k = 2L, tol = 0.5, solver = "highs")
  }
  og1 <- tibble::tibble(team_id = "t1", tract_id = "A")
  expect_message(run(og1), "Computing and caching")
  expect_message(run(og1), "Using cached")
  expect_message(run(tibble::tibble(team_id = "t1", tract_id = "B")), "Computing and caching")
  expect_message(run(NULL), "Computing and caching")
})

test_that("walk_distances reach the roteiro candidates", {
  skip_if_not_installed("ROI.plugin.highs")
  x <- toy_instance()
  # every pair 5 km apart on foot: no roteiro has more than one tract
  pts <- x$points$point_id
  far <- tidyr::expand_grid(origin_id = pts, destination_id = pts) |>
    dplyr::filter(origin_id != destination_id) |>
    dplyr::mutate(distance = 5000)
  res <- surveyzones_assign_teams(x$units, x$points, x$teams, x$eligible, x$durations,
                                  mandatory = x$units$tract_id, k = 3L, tol = 0.5,
                                  walk_distances = far, solver = "highs", use_cache = FALSE)
  expect_true(all(res$roteiros$n_tracts == 1L))
  res_line <- surveyzones_assign_teams(x$units, x$points, x$teams, x$eligible, x$durations,
                                       mandatory = x$units$tract_id, k = 3L, tol = 0.5,
                                       solver = "highs", use_cache = FALSE)
  expect_true(any(res_line$roteiros$n_tracts > 1L))
})

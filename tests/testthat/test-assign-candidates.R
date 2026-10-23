pts_line <- function() {
  # 4 tracts on a parallel, ~300 m apart (0.0028 deg lon at -13 lat ~ 303 m)
  sf::st_as_sf(
    tibble::tibble(
      point_id = c("a1", "a2", "b1", "c1", "c2", "d1"),
      tract_id = c("A", "A", "B", "C", "C", "D"),
      lon = c(-38.5000, -38.5100, -38.4972, -38.4944, -38.4800, -38.4800),
      lat = -13
    ),
    coords = c("lon", "lat"), crs = 4326
  )
}
elig_all <- tidyr::expand_grid(tract_id = c("A", "B", "C", "D"), team_id = c("t1", "t2"))

test_that("one nearby household brings the whole tract", {
  cand <- surveyzones_roteiro_candidates(pts_line(), elig_all, tract_ids = c("A", "B", "C", "D"), radius_m = 500)
  sets <- purrr::map_chr(cand$tract_ids, \(x) paste(x, collapse = "+"))
  expect_true("A+B+C" %in% sets)
  expect_false(any(cand$center_point_id == "b1" & purrr::map_lgl(cand$tract_ids, \(x) "D" %in% x)))
  expect_equal(cand$roteiro_id[1], "R0001")
})

test_that("every tract has a singleton and subsets need a common eligible team", {
  elig <- dplyr::filter(elig_all, !(tract_id == "A" & team_id == "t2"), !(tract_id == "B" & team_id == "t1"))
  cand <- surveyzones_roteiro_candidates(pts_line(), elig, tract_ids = c("A", "B", "C", "D"))
  sets <- purrr::map_chr(cand$tract_ids, \(x) paste(x, collapse = "+"))
  expect_true(all(c("A", "B", "C", "D") %in% sets))
  expect_false("A+B" %in% sets)
  expect_false(anyDuplicated(sets) > 0)
})

test_that("tract without coordinates becomes a warned singleton", {
  cand <- surveyzones_roteiro_candidates(pts_line(), elig_all, tract_ids = c("A", "B", "C", "D", "E"))
  e <- cand[purrr::map_lgl(cand$tract_ids, \(x) identical(x, "E")), ]
  expect_equal(nrow(e), 1L)
  expect_true(is.na(e$center_point_id))
  expect_match(attr(cand, "warnings"), "E", all = FALSE)
})

test_that("max_tracts caps subset size", {
  cand <- surveyzones_roteiro_candidates(pts_line(), elig_all, tract_ids = c("A", "B", "C", "D"), max_tracts = 2L)
  expect_true(all(lengths(cand$tract_ids) <= 2L))
})

test_that("invalid inputs abort", {
  expect_error(surveyzones_roteiro_candidates(pts_line(), elig_all, c("A", "B"), radius_m = 0), "radius_m")
  expect_error(
    surveyzones_roteiro_candidates(dplyr::select(pts_line(), -point_id), elig_all, c("A", "B")),
    "point_id"
  )
})

test_that("teams, span_m, ordering and ids are right", {
  cand <- surveyzones_roteiro_candidates(pts_line(), elig_all, tract_ids = c("A", "B", "C", "D"))
  key <- purrr::map_chr(cand$tract_ids, \(x) paste(x, collapse = "+"))
  abc <- cand[key == "A+B+C", ]
  expect_equal(abc$center_point_id, "b1")
  expect_equal(abc$teams[[1]], c("t1", "t2"))
  expect_equal(abc$span_m, 303, tolerance = 5 / 303)
  expect_true(all(cand$span_m[lengths(cand$tract_ids) == 1L] == 0))
  expect_equal(cand$roteiro_id, sprintf("R%04d", seq_len(nrow(cand))))
  expect_false(is.unsorted(cand$center_tract_id))
})

test_that("dedup keeps smallest span, ties by point_id", {
  pts <- sf::st_as_sf(
    tibble::tibble(
      point_id = c("p2", "p1", "q1"), tract_id = c("A", "A", "B"),
      lon = c(-38.5, -38.5, -38.4972), lat = -13
    ),
    coords = c("lon", "lat"), crs = 4326
  )
  el <- tidyr::expand_grid(tract_id = c("A", "B"), team_id = "t1")
  cand <- surveyzones_roteiro_candidates(pts, el, c("A", "B"))
  ab <- cand[purrr::map_lgl(cand$tract_ids, \(x) identical(x, c("A", "B"))), ]
  expect_equal(nrow(ab), 1L)
  expect_equal(ab$center_point_id, "p1")
  expect_equal(cand$center_point_id[purrr::map_lgl(cand$tract_ids, \(x) identical(x, "A"))], "p1")
})

test_that("tract without eligible team gets singleton with empty teams and warning", {
  el <- dplyr::filter(elig_all, tract_id != "D")
  cand <- surveyzones_roteiro_candidates(pts_line(), el, c("A", "B", "C", "D"))
  d <- cand[purrr::map_lgl(cand$tract_ids, \(x) identical(x, "D")), ]
  expect_equal(nrow(d), 1L)
  expect_equal(d$teams[[1]], character())
  expect_match(attr(cand, "warnings"), "tract D has no eligible team", all = FALSE)
  expect_false(any(purrr::map_lgl(cand$tract_ids, \(x) length(x) > 1L && "D" %in% x)))
})

test_that("bad points abort", {
  expect_error(
    surveyzones_roteiro_candidates(pts_line(), elig_all, c("A", "B", "C")),
    "D"
  )
  p <- pts_line()
  p$tract_id[1] <- NA
  expect_error(surveyzones_roteiro_candidates(p, elig_all, c("A", "B", "C", "D")), "NA")
  p <- pts_line()
  p$point_id[2] <- "a1"
  expect_error(surveyzones_roteiro_candidates(p, elig_all, c("A", "B", "C", "D")), "point_id")
  expect_error(
    surveyzones_roteiro_candidates(pts_line(), elig_all, c("A", "B", "C", "D"), max_tracts = 2.5),
    "max_tracts"
  )
})

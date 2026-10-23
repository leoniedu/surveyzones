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

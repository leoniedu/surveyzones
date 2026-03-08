# Min-Pairwise Zone Distance Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Replace center-to-center zone distances with minimum pairwise tract distances in `surveyzones_sequence_zones()`, so that adjacent zones get low inter-zone distances and better sequence ordering.

**Architecture:** Add a helper `.build_zone_dist_matrix()` that computes `d(zone_A, zone_B) = min(d(tract_i, tract_j))` from the sparse distance table joined with zone assignments. Refactor `surveyzones_sequence_zones()` to build a zone distance matrix directly and call `seriation::seriate()` on it, instead of routing through `.seriate_1d()` with center IDs. Haversine fallback fills missing zone pairs using all access points per zone.

**Tech Stack:** R, dplyr, seriation, sf (for haversine fallback)

---

### Task 1: Write failing test for `.build_zone_dist_matrix()`

**Files:**
- Modify: `tests/testthat/test-sequencing.R`

**Step 1: Write the failing test**

Add this test at the end of `tests/testthat/test-sequencing.R`:

```r
test_that(".build_zone_dist_matrix computes min pairwise distances", {
  # 4 tracts in 2 zones: Z1 = {T1, T2}, Z2 = {T3, T4}
  # T2 and T3 are close (distance 1), T1-T4 are far
  assignments <- tibble::tibble(
    tract_id = paste0("T", 1:4),
    zone_id = c("Z1", "Z1", "Z2", "Z2"),
    partition_id = "P1"
  )

  sparse_dists <- tibble::tibble(
    origin_id = c("T1", "T2", "T1", "T3", "T2", "T4", "T3", "T4",
                  "T2", "T1", "T3", "T1", "T4", "T2", "T4", "T3"),
    destination_id = c("T2", "T1", "T3", "T1", "T3", "T2", "T4", "T3",
                       "T4", "T4", "T2", "T2", "T1", "T3", "T1", "T4"),
    distance = c(2, 2, 10, 10, 1, 8, 3, 3,
                 8, 9, 1, 2, 9, 3, 10, 3)
  )
  # Cross-zone pairs: T1-T3=10, T2-T3=1, T1-T4=9, T2-T4=8
  # min(Z1, Z2) should be 1 (from T2-T3)

  mat <- surveyzones:::.build_zone_dist_matrix(
    zone_ids = c("Z1", "Z2"),
    assignments = assignments,
    sparse_distances = sparse_dists
  )

  expect_equal(nrow(mat), 2)
  expect_equal(ncol(mat), 2)
  expect_equal(rownames(mat), c("Z1", "Z2"))
  expect_equal(mat["Z1", "Z2"], 1)
  expect_equal(mat["Z2", "Z1"], 1)
  expect_equal(mat["Z1", "Z1"], 0)
})
```

**Step 2: Run test to verify it fails**

Run: `Rscript -e 'devtools::test(filter = "sequencing")'`
Expected: FAIL — `.build_zone_dist_matrix` not found.

---

### Task 2: Implement `.build_zone_dist_matrix()`

**Files:**
- Modify: `R/sequencing.R`

**Step 1: Add the helper function**

Insert before `surveyzones_sequence_zones()` in `R/sequencing.R`:

```r
#' Build a Zone-to-Zone Distance Matrix Using Min Pairwise Distances
#'
#' Computes `d(zone_A, zone_B) = min(d(tract_i, tract_j))` for all
#' `i in zone_A, j in zone_B` using the sparse distance table.
#'
#' @param zone_ids Character vector of zone IDs.
#' @param assignments Tibble with `tract_id` and `zone_id` columns.
#' @param sparse_distances Sparse distance tibble.
#'
#' @return A symmetric numeric matrix with zone_ids as row/col names.
#'   Diagonal is 0. Missing zone pairs get `Inf`.
#' @keywords internal
.build_zone_dist_matrix <- function(zone_ids, assignments, sparse_distances) {
  n <- length(zone_ids)
  mat <- matrix(Inf, nrow = n, ncol = n, dimnames = list(zone_ids, zone_ids))
  diag(mat) <- 0

  if (n <= 1) return(mat)

  # Map tract_id -> zone_id
  tract_zone <- assignments |>
    dplyr::select("tract_id", "zone_id")

  # Join zone_ids onto sparse distances
  cross <- sparse_distances |>
    dplyr::inner_join(tract_zone, by = c("origin_id" = "tract_id")) |>
    dplyr::rename(zone_origin = "zone_id") |>
    dplyr::inner_join(tract_zone, by = c("destination_id" = "tract_id")) |>
    dplyr::rename(zone_dest = "zone_id") |>
    dplyr::filter(
      zone_origin != zone_dest,
      zone_origin %in% zone_ids,
      zone_dest %in% zone_ids
    )

  if (nrow(cross) == 0) return(mat)

  # Min distance per zone pair
  zone_dists <- cross |>
    dplyr::summarise(
      min_dist = min(distance),
      .by = c(zone_origin, zone_dest)
    )

  # Fill matrix
  mat[cbind(zone_dists$zone_origin, zone_dists$zone_dest)] <- zone_dists$min_dist

  # Symmetrise: take min of d(A,B) and d(B,A)
  mat <- pmin(mat, t(mat))

  mat
}
```

**Step 2: Run test to verify it passes**

Run: `Rscript -e 'devtools::test(filter = "sequencing")'`
Expected: All tests PASS including the new one.

**Step 3: Commit**

```bash
git add R/sequencing.R tests/testthat/test-sequencing.R
git commit -m "feat: add .build_zone_dist_matrix() for min pairwise zone distances"
```

---

### Task 3: Write failing test for haversine fallback on missing zone pairs

**Files:**
- Modify: `tests/testthat/test-sequencing.R`

**Step 1: Write the failing test**

```r
test_that(".complete_zone_dist_matrix fills Inf pairs via haversine", {
  # 3 zones, but only Z1-Z2 have cross-zone tract distances
  mat <- matrix(Inf, 3, 3, dimnames = list(paste0("Z", 1:3), paste0("Z", 1:3)))
  diag(mat) <- 0
  mat["Z1", "Z2"] <- 5
  mat["Z2", "Z1"] <- 5
  # Z1-Z3 and Z2-Z3 are Inf (missing)

  assignments <- tibble::tibble(
    tract_id = paste0("T", 1:6),
    zone_id = rep(paste0("Z", 1:3), each = 2),
    partition_id = "P1"
  )

  pts <- sf::st_as_sf(
    data.frame(
      tract_id = paste0("T", 1:6),
      lon = c(-38.50, -38.501, -38.502, -38.503, -38.510, -38.511),
      lat = c(-13.00, -13.001, -13.002, -13.003, -13.010, -13.011)
    ),
    coords = c("lon", "lat"), crs = 4326
  )

  result <- surveyzones:::.complete_zone_dist_matrix(
    mat, assignments, pts, speed_kmh = 0.1
  )

  # Z1-Z2 should still be 5 (not overwritten)
  expect_equal(result["Z1", "Z2"], 5)
  # Z1-Z3 and Z2-Z3 should now be finite
  expect_true(is.finite(result["Z1", "Z3"]))
  expect_true(is.finite(result["Z2", "Z3"]))
  # Should be symmetric
  expect_equal(result["Z1", "Z3"], result["Z3", "Z1"])
})
```

**Step 2: Run to verify it fails**

Run: `Rscript -e 'devtools::test(filter = "sequencing")'`
Expected: FAIL — `.complete_zone_dist_matrix` not found.

---

### Task 4: Implement `.complete_zone_dist_matrix()`

**Files:**
- Modify: `R/sequencing.R`

**Step 1: Add the helper function**

Insert after `.build_zone_dist_matrix()`:

```r
#' Fill Missing Zone Pairs via Haversine
#'
#' For zone pairs where `mat[i,j] == Inf`, computes the minimum haversine
#' distance between any access point in zone i and any access point in zone j,
#' converted to travel time using `speed_kmh`.
#'
#' @param mat Numeric matrix from [.build_zone_dist_matrix()].
#' @param assignments Tibble with `tract_id` and `zone_id`.
#' @param access_points sf object with `tract_id` column and POINT geometry.
#' @param speed_kmh Travel speed for haversine conversion (km/h).
#'
#' @return The same matrix with Inf values replaced by haversine-based times.
#' @keywords internal
.complete_zone_dist_matrix <- function(mat, assignments, access_points,
                                       speed_kmh = 0.1) {
  zone_ids <- rownames(mat)
  n <- length(zone_ids)

  # Find missing pairs (upper triangle only, then symmetrise)
  missing_pairs <- which(mat == Inf & upper.tri(mat), arr.ind = TRUE)
  if (nrow(missing_pairs) == 0) return(mat)

  # Map tract -> zone
  tract_zone <- assignments |>
    dplyr::filter(.data$zone_id %in% zone_ids) |>
    dplyr::select("tract_id", "zone_id")

  # Join access_points with zone info
  pts_with_zone <- access_points |>
    dplyr::inner_join(tract_zone, by = "tract_id")

  n_filled <- 0L
  for (k in seq_len(nrow(missing_pairs))) {
    zi <- zone_ids[missing_pairs[k, 1]]
    zj <- zone_ids[missing_pairs[k, 2]]

    pts_i <- pts_with_zone |> dplyr::filter(.data$zone_id == zi)
    pts_j <- pts_with_zone |> dplyr::filter(.data$zone_id == zj)

    if (nrow(pts_i) == 0 || nrow(pts_j) == 0) next

    # Cross-distance matrix (meters), take minimum
    d_mat <- sf::st_distance(pts_i, pts_j)
    min_km <- as.numeric(min(d_mat)) / 1000
    min_time <- min_km / speed_kmh * 60

    mat[zi, zj] <- min_time
    mat[zj, zi] <- min_time
    n_filled <- n_filled + 1L
  }

  if (n_filled > 0) {
    cli::cli_alert_info(
      "Filled {n_filled} missing zone pair{?s} using haversine (speed = {speed_kmh} km/h)."
    )
  }

  mat
}
```

**Step 2: Run tests**

Run: `Rscript -e 'devtools::test(filter = "sequencing")'`
Expected: All PASS.

**Step 3: Commit**

```bash
git add R/sequencing.R tests/testthat/test-sequencing.R
git commit -m "feat: add .complete_zone_dist_matrix() for haversine fallback"
```

---

### Task 5: Write failing test for updated `surveyzones_sequence_zones()`

**Files:**
- Modify: `tests/testthat/test-sequencing.R`

**Step 1: Write the failing test**

This test has two zones where centers are far apart but border tracts are adjacent.
With center-to-center distances, the ordering might be wrong. With min-pairwise, it should be correct.

```r
test_that("sequence_zones uses min pairwise distances, not center-to-center", {
  # Z1 center at C1 (far left), Z2 center at C4 (far right)
  # But Z1 contains T2 (near right) and Z2 contains T3 (near left)
  # So min-pairwise(Z1,Z2) = d(T2,T3) which is small
  # Z3 center at C5 (middle), but all its tracts are in the middle
  #
  # Layout: C1...T2|T3...C4   C5
  # Center-to-center: d(Z1,Z3) < d(Z1,Z2) because C5 is closer to C1 than C4
  # Min-pairwise: d(Z1,Z2) < d(Z1,Z3) because T2-T3 are touching

  assignments <- tibble::tibble(
    tract_id = c("C1", "T2", "T3", "C4", "C5", "T6"),
    zone_id = c("Z1", "Z1", "Z2", "Z2", "Z3", "Z3"),
    partition_id = "P1"
  )

  zones_tbl <- tibble::tibble(
    zone_id = c("Z1", "Z2", "Z3"),
    partition_id = "P1",
    center_tract_id = c("C1", "C4", "C5"),
    total_workload = 2,
    diameter = 1,
    n_tracts = 2L
  )

  # All 6 tracts along a line: C1(0), T2(4), T3(5), C4(9), C5(15), T6(16)
  positions <- c(C1 = 0, T2 = 4, T3 = 5, C4 = 9, C5 = 15, T6 = 16)
  ids <- names(positions)
  pairs <- tidyr::expand_grid(origin_id = ids, destination_id = ids) |>
    dplyr::filter(origin_id != destination_id) |>
    dplyr::mutate(
      distance = abs(positions[origin_id] - positions[destination_id])
    )

  pts <- sf::st_as_sf(
    data.frame(
      tract_id = ids,
      lon = -38.50 - (positions / 1000),
      lat = rep(-13.00, length(ids))
    ),
    coords = c("lon", "lat"), crs = 4326
  )

  plan <- structure(
    list(
      zones = zones_tbl,
      assignments = assignments,
      zone_sequence = NULL, sequence = NULL,
      parameters = list(), diagnostics = list(),
      access_points = pts
    ),
    class = "surveyzones_plan"
  )

  plan <- surveyzones_sequence_zones(
    plan, pairs, pts,
    method = "SPIN_NH", complete_distances = FALSE
  )

  # Z1 and Z2 should be adjacent in the ordering (min pairwise = 1)
  # Z3 should be at one end (far from both)
  ordered <- plan$zone_sequence |>
    dplyr::arrange(zone_order)

  z1_order <- ordered$zone_order[ordered$zone_id == "Z1"]
  z2_order <- ordered$zone_order[ordered$zone_id == "Z2"]

  # Z1 and Z2 must be consecutive

  expect_equal(abs(z1_order - z2_order), 1)
})
```

**Step 2: Run to verify it fails**

Run: `Rscript -e 'devtools::test(filter = "sequencing")'`
Expected: FAIL — current center-to-center approach puts Z3 between Z1 and Z2.

---

### Task 6: Refactor `surveyzones_sequence_zones()` to use min-pairwise distances

**Files:**
- Modify: `R/sequencing.R`

**Step 1: Rewrite `surveyzones_sequence_zones()`**

Replace the body of `surveyzones_sequence_zones()` (lines 167-239 of current `R/sequencing.R`). Keep the signature and roxygen identical except update the `@param access_points` and `@param complete_distances` docs to mention zone pairs instead of center pairs.

The key changes:
1. Build zone distance matrix via `.build_zone_dist_matrix()` instead of using center IDs with `.seriate_1d()`
2. Fill missing zone pairs via `.complete_zone_dist_matrix()` instead of `surveyzones_complete_distances()`
3. Call `seriation::seriate()` directly on the zone distance matrix
4. Extract ordering from the seriation permutation

```r
surveyzones_sequence_zones <- function(
  plan,
  sparse_distances,
  access_points = NULL,
  speed_kmh = 0.1,
  complete_distances = TRUE,
  method = "SPIN_NH",
  control = NULL,
  by_partition = TRUE
) {
  if (!inherits(plan, "surveyzones_plan")) {
    cli::cli_abort("{.arg plan} must be a {.cls surveyzones_plan} object.")
  }

  access_points <- access_points %||% plan$access_points

  # Apply SPIN_NH defaults
  if (method == "SPIN_NH" && is.null(control)) {
    control <- list(sigma = 1)
  }

  # Determine grouping column
  group_col <- if (by_partition) "partition_id" else ".group"
  zones <- plan$zones
  if (!by_partition) {
    zones <- dplyr::mutate(zones, .group = "all")
  }

  zone_lookup <- plan$zones |>
    dplyr::select("zone_id", "partition_id")

  plan$zone_sequence <- zones |>
    tidyr::nest(data = -dplyr::all_of(group_col)) |>
    dplyr::mutate(
      ordered = purrr::map(data, \(z) {
        zone_ids <- z$zone_id
        n <- length(zone_ids)
        if (n <= 1) return(zone_ids)

        # Build min-pairwise zone distance matrix
        mat <- .build_zone_dist_matrix(
          zone_ids, plan$assignments, sparse_distances
        )

        # Fill missing zone pairs via haversine
        if (complete_distances && any(is.infinite(mat[mat != 0]))) {
          if (is.null(access_points)) {
            cli::cli_abort(c(
              "{.arg access_points} is required when {.arg complete_distances} is {.val TRUE}.",
              "i" = "Either pass it explicitly, store it in the plan via {.fun surveyzones_build_zones}, or set {.code complete_distances = FALSE}."
            ))
          }
          mat <- .complete_zone_dist_matrix(
            mat, plan$assignments, access_points, speed_kmh = speed_kmh
          )
        }

        # Seriate
        o <- seriation::seriate(stats::as.dist(mat), method = method,
                                control = control)
        perm <- seriation::get_order(o)
        zone_ids[perm]
      })
    ) |>
    dplyr::select(-data) |>
    tidyr::unnest_longer(ordered, values_to = "zone_id") |>
    dplyr::group_by(dplyr::across(dplyr::all_of(group_col))) |>
    dplyr::mutate(zone_order = dplyr::row_number()) |>
    dplyr::ungroup()

  # When by_partition = FALSE, restore partition_id from the zone lookup
  if (!by_partition) {
    plan$zone_sequence <- plan$zone_sequence |>
      dplyr::left_join(zone_lookup, by = "zone_id") |>
      dplyr::select("partition_id", "zone_id", "zone_order")
  } else {
    plan$zone_sequence <- plan$zone_sequence |>
      dplyr::select("partition_id", "zone_id", "zone_order")
  }

  plan
}
```

**Step 2: Update roxygen for the changed parameters**

Update these `@param` lines:
- `access_points`: "Used to fill missing pairwise distances between zones via haversine" (not "between zone centers")
- `complete_distances`: "missing pairwise distances between zones are filled" (not "between zone centers")

**Step 3: Run all tests**

Run: `Rscript -e 'devtools::test(filter = "sequencing")'`
Expected: All PASS.

**Step 4: Commit**

```bash
git add R/sequencing.R tests/testthat/test-sequencing.R
git commit -m "feat: use min pairwise tract distances for zone sequencing"
```

---

### Task 7: Update existing tests that assume center-to-center behavior

**Files:**
- Modify: `tests/testthat/test-sequencing.R`

**Step 1: Review and fix existing tests**

The existing tests at lines 73-131 ("Spectral seriation produces correct zone_order") use 1-tract-per-zone setups where center == only tract. These should still pass because min-pairwise degenerates to center-to-center when each zone has one tract.

The test at lines 241-288 ("zone sequencing completes incomplete distances via haversine") uses 1-tract-per-zone with incomplete distances. This test's `complete_distances` flow changed — it now uses `.complete_zone_dist_matrix()` instead of `surveyzones_complete_distances()`. Verify it passes; if not, adjust the test data so that the zone distance matrix has `Inf` entries that trigger the haversine fallback.

**Step 2: Run full test suite**

Run: `Rscript -e 'devtools::test(filter = "sequencing")'`
Expected: All PASS.

**Step 3: Run devtools::document() and full check**

Run: `Rscript -e 'devtools::document()'`
Run: `Rscript -e 'devtools::test()'`

**Step 4: Commit if any test adjustments were needed**

```bash
git add tests/testthat/test-sequencing.R
git commit -m "test: update zone sequencing tests for min-pairwise distances"
```

---

### Task 8: Update vignette references

**Files:**
- Modify: `vignettes/articles/census-example.Rmd`

**Step 1: Check if any text references "center-to-center" or "zone centers" in the sequencing sections**

Search for mentions of center distances in the sequencing discussion and update to reflect that zone ordering now uses minimum pairwise tract distances.

**Step 2: Commit if changes were needed**

```bash
git add vignettes/articles/census-example.Rmd
git commit -m "docs: update vignette for min-pairwise zone distances"
```

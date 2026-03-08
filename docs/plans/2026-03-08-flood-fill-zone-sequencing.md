# Flood-Fill Zone Sequencing Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Replace seriation-based zone sequencing with a two-level flood-fill + TSP approach that groups nearby zones by a distance threshold and orders groups via TSP.

**Architecture:** A new internal helper `.flood_fill_groups()` partitions zones into threshold-connected components. `surveyzones_sequence_zones()` is rewritten to: (1) build the zone distance matrix, (2) flood-fill into groups, (3) TSP-seriate within each group, (4) TSP-seriate the groups, (5) concatenate into final zone order. The `method`/`control` params are removed; `threshold` is added.

**Tech Stack:** R, dplyr, seriation (for within/between-group TSP via `.seriate_1d()`), sf

---

### Task 1: Write test and implement `.flood_fill_groups()`

**Files:**
- Modify: `tests/testthat/test-sequencing.R`
- Modify: `R/sequencing.R`

**Step 1: Write the test**

Append to `tests/testthat/test-sequencing.R`:

```r
test_that(".flood_fill_groups finds threshold-connected components", {
  # 5 zones: Z1-Z2 close (d=2), Z2-Z3 close (d=3), Z4-Z5 close (d=1)
  # Z3-Z4 far (d=20)
  # threshold=5 → two groups: {Z1,Z2,Z3} and {Z4,Z5}
  mat <- matrix(Inf, 5, 5, dimnames = list(paste0("Z", 1:5), paste0("Z", 1:5)))
  diag(mat) <- 0
  mat["Z1", "Z2"] <- 2; mat["Z2", "Z1"] <- 2
  mat["Z2", "Z3"] <- 3; mat["Z3", "Z2"] <- 3
  mat["Z1", "Z3"] <- 5; mat["Z3", "Z1"] <- 5
  mat["Z3", "Z4"] <- 20; mat["Z4", "Z3"] <- 20
  mat["Z4", "Z5"] <- 1; mat["Z5", "Z4"] <- 1
  mat["Z1", "Z4"] <- 22; mat["Z4", "Z1"] <- 22
  mat["Z1", "Z5"] <- 23; mat["Z5", "Z1"] <- 23
  mat["Z2", "Z4"] <- 20; mat["Z4", "Z2"] <- 20
  mat["Z2", "Z5"] <- 21; mat["Z5", "Z2"] <- 21
  mat["Z3", "Z5"] <- 21; mat["Z5", "Z3"] <- 21

  groups <- surveyzones:::.flood_fill_groups(mat, threshold = 5)

  expect_equal(length(groups), 2)
  # Each group is a character vector of zone IDs
  group_sizes <- sort(vapply(groups, length, integer(1)))
  expect_equal(group_sizes, c(2L, 3L))

  # Find which group has 3 members and which has 2
  big <- groups[[which(vapply(groups, length, integer(1)) == 3)]]
  small <- groups[[which(vapply(groups, length, integer(1)) == 2)]]
  expect_equal(sort(big), c("Z1", "Z2", "Z3"))
  expect_equal(sort(small), c("Z4", "Z5"))
})

test_that(".flood_fill_groups puts all zones in one group when threshold is large", {
  mat <- matrix(5, 3, 3, dimnames = list(paste0("Z", 1:3), paste0("Z", 1:3)))
  diag(mat) <- 0

  groups <- surveyzones:::.flood_fill_groups(mat, threshold = 10)
  expect_equal(length(groups), 1)
  expect_equal(sort(groups[[1]]), paste0("Z", 1:3))
})

test_that(".flood_fill_groups makes singleton groups when threshold is tiny", {
  mat <- matrix(5, 3, 3, dimnames = list(paste0("Z", 1:3), paste0("Z", 1:3)))
  diag(mat) <- 0

  groups <- surveyzones:::.flood_fill_groups(mat, threshold = 1)
  expect_equal(length(groups), 3)
  all_zones <- sort(unlist(groups))
  expect_equal(all_zones, paste0("Z", 1:3))
})
```

**Step 2: Implement `.flood_fill_groups()`**

Insert in `R/sequencing.R` after `.build_zone_dist_matrix()` and before
the `surveyzones_sequence_zones()` roxygen block:

```r
#' Flood-Fill Zones into Threshold-Connected Groups
#'
#' Partitions zones into groups where every zone in a group is reachable
#' from every other zone in the same group via a chain of pairwise
#' distances <= `threshold`.
#'
#' @param mat Symmetric numeric distance matrix with zone IDs as row/col names.
#' @param threshold Numeric scalar.  Maximum distance for two zones to be
#'   considered connected.
#'
#' @return A list of character vectors, each containing the zone IDs in one
#'   group.
#' @keywords internal
.flood_fill_groups <- function(mat, threshold) {
  zone_ids <- rownames(mat)
  n <- length(zone_ids)
  visited <- rep(FALSE, n)
  names(visited) <- zone_ids
  groups <- list()

  for (seed in zone_ids) {
    if (visited[seed]) next

    # BFS from seed
    queue <- seed
    group <- character(0)
    while (length(queue) > 0) {
      current <- queue[1]
      queue <- queue[-1]
      if (visited[current]) next
      visited[current] <- TRUE
      group <- c(group, current)

      # Find unvisited neighbors within threshold
      neighbors <- zone_ids[!visited & mat[current, ] <= threshold & mat[current, ] > 0]
      queue <- c(queue, neighbors)
    }

    groups <- c(groups, list(group))
  }

  groups
}
```

**Step 3: Run tests**

Run: `Rscript -e 'devtools::test(filter = "sequencing")'`
Expected: All PASS.

**Step 4: Commit**

```bash
git add R/sequencing.R tests/testthat/test-sequencing.R
git commit -m "feat: add .flood_fill_groups() for threshold-connected zone grouping"
```

---

### Task 2: Rewrite `surveyzones_sequence_zones()` and update tests

**Files:**
- Modify: `R/sequencing.R`
- Modify: `tests/testthat/test-sequencing.R`

**Step 1: Replace `surveyzones_sequence_zones()`**

Replace the entire function (from `#' Sequence Zones Within Each Partition`
roxygen through closing `}`) with:

```r
#' Sequence Zones Within Each Partition
#'
#' Given a solved plan, finds a spatial ordering for the zones within each
#' partition using a flood-fill + TSP approach.  Zones within `threshold`
#' distance of each other are grouped together, then groups and zones within
#' groups are ordered via TSP seriation on minimum pairwise tract distances.
#'
#' @param plan A `surveyzones_plan` object.
#' @param sparse_distances Sparse distance table (as returned by
#'   [surveyzones_compute_sparse_distances()]).
#' @param access_points An sf object with POINT geometries and a
#'   `tract_id` column.  Used to fill missing pairwise distances via
#'   haversine ([surveyzones_complete_distances()]) so that every
#'   cross-zone tract pair has a distance.
#'   Defaults to `plan$access_points` when `NULL`.  Required when
#'   `complete_distances = TRUE`.
#' @param speed_kmh Numeric scalar.  Assumed travel speed for haversine
#'   fill-in (km/h).  Default `0.1` — intentionally harsh because missing
#'   pairs likely represent real barriers (rivers, mountains).
#' @param complete_distances Logical scalar.  When `TRUE` (default),
#'   missing pairwise distances are filled using haversine via
#'   [surveyzones_complete_distances()].  Set to `FALSE` to skip
#'   imputation (e.g., when distances are already complete).
#' @param threshold Numeric scalar.  Maximum distance (in the same units as
#'   `sparse_distances`, typically minutes) for two zones to be grouped
#'   together.  Default `10`.
#' @param by_partition Logical scalar.  When `TRUE` (default), zones are
#'   sequenced independently within each partition.  When `FALSE`, all zones
#'   are sequenced together ignoring partition boundaries.
#'
#' @return The same `plan` with `plan$zone_sequence` populated — a tibble
#'   with columns `partition_id`, `zone_id`, `zone_order`.
#'
#' @export
surveyzones_sequence_zones <- function(
  plan,
  sparse_distances,
  access_points = NULL,
  speed_kmh = 0.1,
  complete_distances = TRUE,
  threshold = 10,
  by_partition = TRUE
) {
  if (!inherits(plan, "surveyzones_plan")) {
    cli::cli_abort("{.arg plan} must be a {.cls surveyzones_plan} object.")
  }

  access_points <- access_points %||% plan$access_points

  # Complete missing distances using haversine so every cross-zone pair exists
  if (complete_distances) {
    if (is.null(access_points)) {
      cli::cli_abort(c(
        "{.arg access_points} is required when {.arg complete_distances} is {.val TRUE}.",
        "i" = "Either pass it explicitly, store it in the plan via {.fun surveyzones_build_zones}, or set {.code complete_distances = FALSE}."
      ))
    }
    sparse_distances <- surveyzones_complete_distances(
      sparse_distances, access_points, speed_kmh = speed_kmh
    )
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

        # Flood-fill into threshold-connected groups
        groups <- .flood_fill_groups(mat, threshold)

        if (length(groups) == 1) {
          # Single group: just TSP-seriate all zones
          return(.seriate_1d(zone_ids, sparse_distances, method = "TSP"))
        }

        # Seriate within each group
        ordered_groups <- purrr::map(groups, \(g) {
          if (length(g) <= 2) return(g)
          .seriate_1d(g, sparse_distances, method = "TSP")
        })

        # Build group-to-group distance matrix
        ng <- length(ordered_groups)
        group_names <- paste0("G", seq_len(ng))
        group_mat <- matrix(Inf, ng, ng, dimnames = list(group_names, group_names))
        diag(group_mat) <- 0

        for (i in seq_len(ng - 1)) {
          for (j in (i + 1):ng) {
            d <- min(mat[groups[[i]], groups[[j]], drop = FALSE])
            group_mat[i, j] <- d
            group_mat[j, i] <- d
          }
        }

        # Seriate groups
        if (ng <= 2) {
          group_order <- seq_len(ng)
        } else {
          group_dist <- stats::as.dist(group_mat)
          o <- seriation::seriate(group_dist, method = "TSP")
          group_order <- seriation::get_order(o)
        }

        # Concatenate: group order × within-group order
        unlist(ordered_groups[group_order])
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

**Step 2: Update existing tests**

The following tests call `surveyzones_sequence_zones()` with `method =`
which no longer exists.  Update them to remove `method`/`control` and
use `threshold` instead.

**Test "Spectral seriation produces correct zone_order" (line ~73):**
Replace the call at line 117:
```r
  # Old:
  plan <- surveyzones_sequence_zones(plan, pairs, pts, method = "Spectral")
  # New:
  plan <- surveyzones_sequence_zones(plan, pairs, pts, threshold = 100)
```
Also update the test name and remove the Spectral-specific assertion
(linear order up to reflection).  The flood-fill with a large threshold
puts all 4 zones in one group, TSP-seriates them, so they should still
be a permutation of Z1-Z4.  Replace the entire test with:
```r
test_that("zone sequencing produces correct zone_order for line topology", {
  # 4 zones in a line: Z1--Z2--Z3--Z4, all within threshold
  zones_tbl <- tibble::tibble(
    zone_id = paste0("Z", 1:4),
    partition_id = "P1",
    center_tract_id = paste0("C", 1:4),
    total_workload = 1,
    diameter = 0,
    n_tracts = 1L
  )

  ids <- paste0("C", 1:4)
  pairs <- tidyr::expand_grid(origin_id = ids, destination_id = ids) |>
    dplyr::filter(origin_id != destination_id) |>
    dplyr::mutate(
      i = as.integer(gsub("C", "", origin_id)),
      j = as.integer(gsub("C", "", destination_id)),
      distance = abs(i - j)
    ) |>
    dplyr::select(origin_id, destination_id, distance)

  pts <- sf::st_as_sf(
    data.frame(
      tract_id = ids,
      lon = -38.50 - (seq_along(ids) - 1) * 0.001,
      lat = -13.00 - (seq_along(ids) - 1) * 0.001
    ),
    coords = c("lon", "lat"), crs = 4326
  )

  plan <- structure(
    list(
      zones = zones_tbl,
      assignments = tibble::tibble(
        tract_id = paste0("C", 1:4),
        zone_id = paste0("Z", 1:4),
        partition_id = "P1"
      ),
      zone_sequence = NULL, sequence = NULL,
      parameters = list(), diagnostics = list()
    ),
    class = "surveyzones_plan"
  )

  plan <- surveyzones_sequence_zones(plan, pairs, pts, threshold = 100)

  expect_false(is.null(plan$zone_sequence))
  expect_true("zone_order" %in% names(plan$zone_sequence))
  expect_equal(nrow(plan$zone_sequence), 4)
  expect_equal(sort(plan$zone_sequence$zone_order), 1:4)

  # TSP on a line should recover the linear order (up to reflection)
  ordered_zones <- plan$zone_sequence |>
    dplyr::arrange(zone_order) |>
    dplyr::pull(zone_id)
  zone_nums <- as.integer(gsub("Z", "", ordered_zones))
  diffs <- diff(zone_nums)
  expect_true(all(diffs > 0) || all(diffs < 0))
})
```

**Test "seriation handles single-zone partition" (line ~133):**
Replace line 166:
```r
  # Old:
  plan <- surveyzones_sequence_zones(plan, pairs, pts, method = "Spectral")
  # New:
  plan <- surveyzones_sequence_zones(plan, pairs, pts)
```

**Test "zone sequencing completes incomplete distances via haversine" (line ~241):**
Replace line 282-284:
```r
  # Old:
  plan <- surveyzones_sequence_zones(
    plan, incomplete_pairs, pts, method = "Spectral"
  )
  # New:
  plan <- surveyzones_sequence_zones(plan, incomplete_pairs, pts)
```

**Test "by_partition = FALSE sequences all zones together" (line ~290):**
Replace line 333-335:
```r
  # Old:
  plan <- surveyzones_sequence_zones(
    plan, pairs, pts, method = "Spectral", by_partition = FALSE
  )
  # New:
  plan <- surveyzones_sequence_zones(
    plan, pairs, pts, threshold = 100, by_partition = FALSE
  )
```

**Test "sequence_zones uses min pairwise distances, not center-to-center" (line ~427):**
Replace line 481-484:
```r
  # Old:
  plan <- surveyzones_sequence_zones(
    plan, pairs, pts,
    method = "SPIN_NH", complete_distances = FALSE
  )
  # New:
  plan <- surveyzones_sequence_zones(
    plan, pairs, pts,
    threshold = 100, complete_distances = FALSE
  )
```

**Step 3: Add a test for flood-fill grouping through the public API**

Append to `tests/testthat/test-sequencing.R`:

```r
test_that("zone sequencing groups zones by threshold", {
  # 6 zones in two clusters: Z1-Z3 close together, Z4-Z6 close together
  # Clusters are far apart (distance 50)
  zones_tbl <- tibble::tibble(
    zone_id = paste0("Z", 1:6),
    partition_id = "P1",
    center_tract_id = paste0("C", 1:6),
    total_workload = 1,
    diameter = 0,
    n_tracts = 1L
  )

  # C1(0), C2(2), C3(4), C4(54), C5(56), C6(58)
  positions <- c(C1 = 0, C2 = 2, C3 = 4, C4 = 54, C5 = 56, C6 = 58)
  ids <- names(positions)
  pairs <- tidyr::expand_grid(origin_id = ids, destination_id = ids) |>
    dplyr::filter(origin_id != destination_id) |>
    dplyr::mutate(
      distance = abs(positions[origin_id] - positions[destination_id])
    )

  pts <- sf::st_as_sf(
    data.frame(
      tract_id = ids,
      lon = -38.50 - (positions / 10000),
      lat = rep(-13.00, length(ids))
    ),
    coords = c("lon", "lat"), crs = 4326
  )

  plan <- structure(
    list(
      zones = zones_tbl,
      assignments = tibble::tibble(
        tract_id = ids,
        zone_id = paste0("Z", 1:6),
        partition_id = "P1"
      ),
      zone_sequence = NULL, sequence = NULL,
      parameters = list(), diagnostics = list()
    ),
    class = "surveyzones_plan"
  )

  plan <- surveyzones_sequence_zones(
    plan, pairs, pts, threshold = 10, complete_distances = FALSE
  )

  ordered <- plan$zone_sequence |> dplyr::arrange(zone_order)

  # Z1-Z3 should be consecutive and Z4-Z6 should be consecutive
  cluster1_orders <- ordered$zone_order[ordered$zone_id %in% paste0("Z", 1:3)]
  cluster2_orders <- ordered$zone_order[ordered$zone_id %in% paste0("Z", 4:6)]

  # Each cluster occupies 3 consecutive positions
  expect_equal(diff(sort(cluster1_orders)), c(1L, 1L))
  expect_equal(diff(sort(cluster2_orders)), c(1L, 1L))

  # The two clusters don't interleave
  expect_true(
    max(cluster1_orders) < min(cluster2_orders) ||
    max(cluster2_orders) < min(cluster1_orders)
  )
})
```

**Step 4: Run all tests**

Run: `Rscript -e 'devtools::test(filter = "sequencing")'`
Expected: All PASS.

**Step 5: Commit**

```bash
git add R/sequencing.R tests/testthat/test-sequencing.R
git commit -m "feat: replace seriation with flood-fill + TSP for zone sequencing"
```

---

### Task 3: Update docs and run full check

**Files:**
- Modify: `vignettes/articles/census-example.Rmd`

**Step 1: Update vignette comment**

Replace line 226:
```r
# Old:
# Zone visit order via SPIN_NH seriation on min pairwise tract distances
# New:
# Zone visit order via flood-fill + TSP on min pairwise tract distances
```

**Step 2: Run devtools::document() and full test suite**

Run: `Rscript -e 'devtools::document()'`
Run: `Rscript -e 'devtools::test()'`

**Step 3: Commit**

```bash
git add R/sequencing.R vignettes/articles/census-example.Rmd man/
git commit -m "docs: update docs for flood-fill zone sequencing"
```

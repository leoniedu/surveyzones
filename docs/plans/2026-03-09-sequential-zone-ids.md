# Sequential Zone IDs Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Replace opaque center-tract zone IDs with sequential `{partition_id}_{group}.{pos:03d}` identifiers applied automatically at the end of `surveyzones_sequence()`.

**Architecture:** A new internal helper `.rename_zones(plan)` builds an old→new lookup from `zone_sequence` and applies it to all four tables (`assignments`, `zones`, `zone_sequence`, `sequence`). It also drops the now-redundant `group_id` column from `assignments`. Called as the last step of `surveyzones_sequence()`. Follow-up: revert the `.travel_group` workaround in `pns.zonas`.

**Tech Stack:** R, dplyr, testthat 3, devtools

---

### Task 1: Write failing tests for `.rename_zones`

**Files:**
- Modify: `tests/testthat/test-sequencing.R` (append at end)

**Step 1: Append the tests**

```r
# ── .rename_zones ──────────────────────────────────────────────────────────────

test_that(".rename_zones produces correct ID format", {
  plan <- structure(
    list(
      zones = tibble::tibble(
        zone_id = c("Z1", "Z2", "Z3"),
        partition_id = "P1",
        center_tract_id = c("C1", "C2", "C3"),
        total_workload = 1, diameter = 0, n_tracts = 1L
      ),
      assignments = tibble::tibble(
        tract_id = c("T1", "T2", "T3"),
        zone_id  = c("Z1", "Z2", "Z3"),
        partition_id = "P1"
      ),
      zone_sequence = tibble::tibble(
        partition_id = "P1",
        zone_id      = c("Z1", "Z2", "Z3"),
        zone_order   = 1:3,
        group_id     = c(1L, 1L, 2L)   # Z1+Z2 in group 1, Z3 in group 2
      ),
      sequence = tibble::tibble(
        zone_id     = c("Z1", "Z2", "Z3"),
        tract_id    = c("T1", "T2", "T3"),
        visit_order = 1L
      ),
      parameters = list(), diagnostics = list()
    ),
    class = "surveyzones_plan"
  )

  out <- surveyzones:::.rename_zones(plan)

  # zone_sequence uses new IDs
  expect_equal(sort(out$zone_sequence$zone_id), c("P1_1.001", "P1_1.002", "P1_2.001"))
  # assignments uses new IDs
  expect_equal(sort(out$assignments$zone_id), c("P1_1.001", "P1_1.002", "P1_2.001"))
  # zones uses new IDs
  expect_equal(sort(out$zones$zone_id), c("P1_1.001", "P1_1.002", "P1_2.001"))
  # sequence uses new IDs
  expect_equal(sort(out$sequence$zone_id), c("P1_1.001", "P1_1.002", "P1_2.001"))
  # group_id dropped from assignments
  expect_false("group_id" %in% names(out$assignments))
  # center_tract_id preserved in zones
  expect_true("center_tract_id" %in% names(out$zones))
})

test_that(".rename_zones works across two partitions", {
  plan <- structure(
    list(
      zones = tibble::tibble(
        zone_id = c("Z1", "Z2", "Z3", "Z4"),
        partition_id = c("PA", "PA", "PB", "PB"),
        center_tract_id = paste0("C", 1:4),
        total_workload = 1, diameter = 0, n_tracts = 1L
      ),
      assignments = tibble::tibble(
        tract_id     = paste0("T", 1:4),
        zone_id      = c("Z1", "Z2", "Z3", "Z4"),
        partition_id = c("PA", "PA", "PB", "PB"),
        group_id     = c("g1", "g1", "g2", "g2")  # char group_id to drop
      ),
      zone_sequence = tibble::tibble(
        partition_id = c("PA", "PA", "PB", "PB"),
        zone_id      = c("Z1", "Z2", "Z3", "Z4"),
        zone_order   = c(1L, 2L, 1L, 2L),
        group_id     = c(1L, 1L, 1L, 1L)
      ),
      sequence = tibble::tibble(
        zone_id     = c("Z1", "Z2", "Z3", "Z4"),
        tract_id    = paste0("T", 1:4),
        visit_order = 1L
      ),
      parameters = list(), diagnostics = list()
    ),
    class = "surveyzones_plan"
  )

  out <- surveyzones:::.rename_zones(plan)

  pa_ids <- sort(out$zone_sequence$zone_id[out$zone_sequence$partition_id == "PA"])
  pb_ids <- sort(out$zone_sequence$zone_id[out$zone_sequence$partition_id == "PB"])
  expect_equal(pa_ids, c("PA_1.001", "PA_1.002"))
  expect_equal(pb_ids, c("PB_1.001", "PB_1.002"))
  # IDs are globally unique
  expect_equal(length(unique(out$zone_sequence$zone_id)), 4L)
  # group_id dropped from assignments
  expect_false("group_id" %in% names(out$assignments))
})
```

**Step 2: Run tests to confirm they fail**

```bash
cd /Users/eleon/gitlab/ibgeba/surveyzones
Rscript -e "devtools::test(filter = 'sequencing')" 2>&1 | grep -E "FAIL|Error|rename"
```

Expected: `Error ... .rename_zones not found` or similar.

---

### Task 2: Implement `.rename_zones`

**Files:**
- Modify: `R/sequencing.R` (append helper before the final `surveyzones_sequence_tracts` function)

**Step 1: Add the helper**

Insert after `.seriate_quietly` and before `surveyzones_sequence_tracts` (around line 228):

```r
#' Rename zone IDs to sequential {partition_id}_{group}.{pos:003d} format
#'
#' Applied at the end of [surveyzones_sequence()] after zone_sequence is
#' computed.  Builds an old→new lookup and applies it to all four plan tables.
#' Also drops `group_id` from `assignments` (character phase-1 group, now
#' redundant with `zones$center_tract_id`).
#'
#' @param plan A `surveyzones_plan` with `zone_sequence` already set.
#' @return The plan with renamed zone_ids in all four tables.
#' @keywords internal
.rename_zones <- function(plan) {
  # Build lookup: pos = rank within (partition_id, group_id) ordered by zone_order
  lookup <- plan$zone_sequence |>
    dplyr::arrange(.data$partition_id, .data$group_id, .data$zone_order) |>
    dplyr::mutate(
      pos = dplyr::row_number(),
      .by = c("partition_id", "group_id")
    ) |>
    dplyr::mutate(
      new_zone_id = paste0(
        .data$partition_id, "_",
        .data$group_id, ".",
        sprintf("%03d", .data$pos)
      )
    ) |>
    dplyr::select("zone_id", "new_zone_id")

  id_map <- stats::setNames(lookup$new_zone_id, lookup$zone_id)

  # Apply to all four tables
  plan$zone_sequence$zone_id <- id_map[plan$zone_sequence$zone_id]
  plan$zones$zone_id         <- id_map[plan$zones$zone_id]
  plan$assignments$zone_id   <- id_map[plan$assignments$zone_id]
  plan$sequence$zone_id      <- id_map[plan$sequence$zone_id]

  # Drop character group_id from assignments (redundant with zones$center_tract_id)
  plan$assignments <- dplyr::select(plan$assignments, -dplyr::any_of("group_id"))

  plan
}
```

**Step 2: Run new tests — expect PASS**

```bash
Rscript -e "devtools::test(filter = 'sequencing')" 2>&1 | grep -E "PASS|FAIL|rename"
```

Expected: both new tests PASS.

**Step 3: Commit**

```bash
git add R/sequencing.R tests/testthat/test-sequencing.R
git commit -m "feat: add .rename_zones helper with tests"
```

---

### Task 3: Wire `.rename_zones` into `surveyzones_sequence()`

**Files:**
- Modify: `R/sequencing.R` lines 77–85

**Step 1: Add the call at end of `surveyzones_sequence()`**

Replace:
```r
  plan$zone_sequence <- .sequence_zones(
    plan, sparse_distances, control, threshold, by_partition
  )

  plan$sequence <- .sequence_tracts_with_orientation(
    plan, sparse_distances, method, control
  )

  plan
```

With:
```r
  plan$zone_sequence <- .sequence_zones(
    plan, sparse_distances, control, threshold, by_partition
  )

  plan$sequence <- .sequence_tracts_with_orientation(
    plan, sparse_distances, method, control
  )

  .rename_zones(plan)
```

**Step 2: Run the full test suite**

```bash
Rscript -e "devtools::test()" 2>&1 | tail -20
```

Expected: some existing sequencing tests FAIL (they reference old zone IDs like `"Z1"`, `"ZA"`). Note which tests fail — fix them in Task 4.

---

### Task 4: Update existing tests broken by the rename

**Files:**
- Modify: `tests/testthat/test-sequencing.R`

The failing tests reference original zone IDs (`"Z1"`, `"ZA"`, `"ZB"`) in `zone_sequence` or `sequence` after calling `surveyzones_sequence()`. After the rename, those IDs no longer exist.

**Step 1: Fix "zone sequencing produces correct zone_order for line topology" (~line 128)**

The test checks that the TSP recovers linear order by parsing zone numbers from IDs with `gsub("Z", ...)`. Replace that check with one that uses only `zone_order`:

```r
  # TSP on a line should recover the linear order (up to reflection)
  # Use zone_order directly — IDs have been renamed and are no longer "Z1"..."Z4"
  ordered_zones <- plan$zone_sequence |>
    dplyr::arrange(zone_order) |>
    dplyr::pull(zone_id)
  # All 4 zones present, in some order
  expect_length(ordered_zones, 4)
  expect_equal(sort(plan$zone_sequence$zone_order), 1:4)
  # All new IDs match the expected format {partition}_{group}.{pos:03d}
  expect_true(all(grepl("^P1_\\d+\\.\\d{3}$", ordered_zones)))
```

**Step 2: Fix "gap-based grouping separates two distant clusters" (~line 557)**

The test identifies cluster zones by original IDs (`paste0("Z", 1:3)`). After rename the IDs are gone; use `group_id` in `zone_sequence` instead:

```r
  ordered <- plan$zone_sequence |> dplyr::arrange(zone_order)

  # Two distinct group_ids should exist (one per cluster)
  groups <- sort(unique(ordered$group_id))
  expect_length(groups, 2)

  cluster1_orders <- ordered$zone_order[ordered$group_id == groups[1]]
  cluster2_orders <- ordered$zone_order[ordered$group_id == groups[2]]

  # Each cluster occupies 3 consecutive positions
  expect_equal(diff(sort(cluster1_orders)), c(1L, 1L))
  expect_equal(diff(sort(cluster2_orders)), c(1L, 1L))

  # The two clusters don't interleave
  expect_true(
    max(cluster1_orders) < min(cluster2_orders) ||
    max(cluster2_orders) < min(cluster1_orders)
  )
```

**Step 3: Fix "tract orientation orients second zone's entry toward first zone's exit" (~line 582)**

This test queries `plan$sequence` by zone_id after sequencing, e.g. `plan$sequence$zone_id == "ZA"`. Replace with a lookup using `zone_sequence` to get the new IDs:

```r
  plan <- surveyzones_sequence(
    plan, pairs, complete_distances = FALSE
  )

  # Look up renamed IDs for the two zones (ZA → first zone, ZB → second)
  # zone_sequence maps old→new; we find them by zone_order
  zseq <- plan$zone_sequence |> dplyr::arrange(zone_order)
  id_first  <- zseq$zone_id[1]
  id_second <- zseq$zone_id[2]

  first_tracts  <- plan$sequence$tract_id[plan$sequence$zone_id == id_first]
  second_tracts <- plan$sequence$tract_id[plan$sequence$zone_id == id_second]
```

*(keep the rest of the orientation assertion logic unchanged — it only uses `first_tracts` and `second_tracts`)*

**Step 4: Fix "sequencing returns a permutation of zone tracts" (~line 37)**

This test loops over `plan$zones$zone_id` and looks up `plan$assignments$zone_id == zid`. Both are renamed consistently, so the join still works — but verify:

```bash
Rscript -e "devtools::test(filter = 'sequencing')" 2>&1 | grep -E "PASS|FAIL"
```

If it still fails, check whether the test manually constructs `plan` with literal zone IDs — if so apply the same pattern as Step 3 (look up renamed IDs from `zone_sequence`).

**Step 5: Run full test suite**

```bash
Rscript -e "devtools::test()" 2>&1 | tail -20
```

Expected: all tests PASS.

**Step 6: Commit**

```bash
git add tests/testthat/test-sequencing.R
git commit -m "test: update sequencing tests for renamed zone IDs"
```

---

### Task 5: Fix `pns.zonas/R/map.R` — revert `.travel_group` workaround

Now that `assignments` no longer has `group_id`, the naming conflict is gone. Simplify back to a direct reference.

**Files:**
- Modify: `/Users/eleon/gitlab/ibgeba/pns.zonas/R/map.R` lines 60–73

**Step 1: Replace the workaround block**

Replace:
```r
  # Color by travel group_id when zone_sequence is available, otherwise by zone_id.
  # Use a distinct rename to avoid conflict: plan$assignments may already have
  # a character group_id (phase-1 group); zone_sequence$group_id is an integer
  # travel group — renaming prevents dplyr from creating group_id.x/.y suffixes.
  if (!is.null(plan$zone_sequence) && "group_id" %in% names(plan$zone_sequence)) {
    group_lookup <- plan$zone_sequence |>
      dplyr::select("zone_id", .travel_group = "group_id") |>
      dplyr::distinct()
    zonas <- dplyr::left_join(zonas, group_lookup, by = "zone_id") |>
      dplyr::mutate(color = color_by(.data$.travel_group)) |>
      dplyr::select(-".travel_group")
  } else {
    zonas <- dplyr::mutate(zonas, color = color_by(.data$zone_id))
  }
```

With:
```r
  # Color by travel group when zone_sequence is available, otherwise by zone_id
  if (!is.null(plan$zone_sequence) && "group_id" %in% names(plan$zone_sequence)) {
    group_lookup <- plan$zone_sequence |>
      dplyr::distinct(.data$zone_id, .data$group_id)
    zonas <- dplyr::left_join(zonas, group_lookup, by = "zone_id") |>
      dplyr::mutate(color = color_by(.data$group_id))
  } else {
    zonas <- dplyr::mutate(zonas, color = color_by(.data$zone_id))
  }
```

**Step 2: Commit**

```bash
cd /Users/eleon/gitlab/ibgeba/pns.zonas
git add R/map.R
git commit -m "fix: simplify group_id coloring now that assignments has no group_id"
```

---

### Task 6: Integration smoke test

**Step 1: Reload both packages and verify end-to-end**

```r
devtools::load_all("/Users/eleon/gitlab/ibgeba/surveyzones")
devtools::load_all("/Users/eleon/gitlab/ibgeba/pns.zonas")

plan <- readr::read_rds("/Users/eleon/gitlab/ibgeba/pns.zonas/vignettes/zonas_test.rds")

# Re-run sequencing on the saved plan (it has old-style IDs)
d <- readr::read_rds("/Users/eleon/gitlab/ibgeba/pns.zonas/vignettes/d_upas_osrm_test.rds")
plan_new <- surveyzones_sequence(plan, d, threshold = 10, complete_distances = FALSE)

# IDs should match new format
head(plan_new$assignments$zone_id)   # e.g. "28_1.001"
head(plan_new$zone_sequence$zone_id) # same
stopifnot(!("group_id" %in% names(plan_new$assignments)))
stopifnot(all(grepl("^\\w+_\\d+\\.\\d{3}$", plan_new$assignments$zone_id)))
cat("OK\n")
```

**Step 2: Run surveyzones full test suite one final time**

```bash
cd /Users/eleon/gitlab/ibgeba/surveyzones
Rscript -e "devtools::test()" 2>&1 | tail -5
```

Expected: all tests pass, 0 failures.

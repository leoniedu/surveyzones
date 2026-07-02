# Zone Grouping and Labeling Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add `surveyzones_group_zones()` — given externally-defined (not solver-built) zones, find each zone's center, group zones within `D_max` by reusing the existing tract-level p-median machinery one level up, and assign stable `ggg.zzz` labels via two independent partition-wide TSPs, plus an oriented tract-level visit sequence.

**Architecture:** New file `R/zone_grouping.R` holds the pipeline: input validation, medoid center-finding, a shared zone-to-zone sparse-distance helper, the reused p-median grouping call, the two-TSP labeling step, and tract-level orientation. A new file `R/zone_groups_plan.R` holds the `surveyzones_zone_groups` S3 constructor and print method, mirroring the existing `surveyzones_plan.R` split. `surveyzones_group_zones()` orchestrates all of this per partition and binds the results together. No existing files are modified.

**Tech Stack:** R, dplyr, tidyr, purrr, tibble, cli, rlang, seriation (via existing `.seriate_1d()`), testthat 3, devtools, roxygen2.

## Global Constraints

- Follow existing code style exactly: base pipe `|>`, dplyr verbs with the `.data` pronoun, `cli::cli_abort()` for errors, `@keywords internal` roxygen tag on internal functions, `@export` on public API, tibbles (not data.frame) for all tabular returns.
- Reuse these existing internal/exported functions **unmodified** — do not edit `R/zones.R`, `R/sequencing.R`, `R/distances.R`, `R/partition.R`: `.build_zone_dist_matrix()`, `.matrix_to_sparse_dt()`, `.seriate_1d()`, `surveyzones_sequence_tracts()`, `surveyzones_build_zones_single()`, `surveyzones_build_zones_single_mem`, `validate_solver()`, `surveyzones_complete_distances()`.
- No new package dependencies — everything needed is already in `DESCRIPTION`'s `Imports`. Do not edit `DESCRIPTION`.
- Design reference: `docs/plans/2026-07-02-zone-grouping-design.md` (includes the `starting_points` / anchoring addendum).
- Starting-point anchoring uses `sf::st_distance()` directly (same pattern as `surveyzones_complete_distances()`'s haversine fill-in) — no new package dependency.
- Every test that triggers an actual solver call must use `skip_if_not_installed("ROI.plugin.glpk")` (mirrors `tests/testthat/test-zones.R`). Trivial single-zone/single-group paths that never reach the solver don't need this skip.
- Always pass `use_cache = FALSE` in tests to avoid polluting the on-disk memoise cache (mirrors `tests/testthat/test-zones-spopt.R`).
- **Environment prerequisite:** this repo's R package dependencies (dplyr, tibble, purrr, tidyr, cli, rlang, sf, seriation, igraph, memoise, cachem, ompr, ompr.roi, ROI, ROI.plugin.glpk, testthat, devtools) must be installed before running any step below. If `Rscript -e 'requireNamespace("devtools", quietly = TRUE)'` prints `FALSE`, install them first — this is one-time environment setup, not part of this plan's task cycle.

---

### Task 1: Input validation — `.validate_zone_assignments()`

**Files:**
- Create: `R/zone_grouping.R`
- Create: `tests/testthat/test-zone-grouping.R`

**Interfaces:**
- Produces: `.validate_zone_assignments(zone_assignments, call = rlang::caller_env())` → tibble with `tract_id` (character), `zone_id` (character), `partition_id` (character, defaults to `"all"` when the input has no `partition_id` column). Aborts on: missing `tract_id`/`zone_id` columns, a `tract_id` appearing under more than one `zone_id`, or a `zone_id` whose rows span more than one `partition_id`.

- [ ] **Step 1: Create the test file with failing tests**

```r
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
```

- [ ] **Step 2: Run tests to confirm they fail**

```bash
Rscript -e "devtools::load_all(quiet = TRUE); testthat::test_file('tests/testthat/test-zone-grouping.R')"
```

Expected: errors like `could not find function ".validate_zone_assignments"` for every test — the file `R/zone_grouping.R` doesn't exist yet.

- [ ] **Step 3: Create `R/zone_grouping.R` with the validation function**

```r
#' Validate and Normalize External Zone Assignments
#'
#' @param zone_assignments A data.frame with `tract_id`, `zone_id`, and
#'   optionally `partition_id`.
#' @param call Caller environment for error reporting.
#' @return A tibble with `tract_id`, `zone_id`, `partition_id` (character),
#'   defaulting `partition_id` to `"all"` when absent from the input.
#' @keywords internal
.validate_zone_assignments <- function(zone_assignments, call = rlang::caller_env()) {
  if (!is.data.frame(zone_assignments)) {
    cli::cli_abort("{.arg zone_assignments} must be a data.frame.", call = call)
  }

  required <- c("tract_id", "zone_id")
  missing_cols <- setdiff(required, names(zone_assignments))
  if (length(missing_cols) > 0) {
    cli::cli_abort(
      "{.arg zone_assignments} is missing column{?s}: {.val {missing_cols}}.",
      call = call
    )
  }

  out <- tibble::tibble(
    tract_id = as.character(zone_assignments$tract_id),
    zone_id = as.character(zone_assignments$zone_id),
    partition_id = if ("partition_id" %in% names(zone_assignments)) {
      as.character(zone_assignments$partition_id)
    } else {
      "all"
    }
  )

  if (anyDuplicated(out$tract_id)) {
    dupes <- unique(out$tract_id[duplicated(out$tract_id)])
    cli::cli_abort(
      "{.field tract_id} values must be unique across zones: {.val {dupes}}.",
      call = call
    )
  }

  multi_partition <- out |>
    dplyr::distinct(.data$zone_id, .data$partition_id) |>
    dplyr::count(.data$zone_id) |>
    dplyr::filter(.data$n > 1) |>
    dplyr::pull(.data$zone_id)

  if (length(multi_partition) > 0) {
    cli::cli_abort(
      "Zone{?s} span more than one partition_id: {.val {multi_partition}}.",
      call = call
    )
  }

  out
}
```

- [ ] **Step 4: Run tests — expect PASS**

```bash
Rscript -e "devtools::load_all(quiet = TRUE); testthat::test_file('tests/testthat/test-zone-grouping.R')"
```

Expected: `[ FAIL 0 | ... ]` — 0 failures. (testthat's PASS count reflects individual `expect_*()` calls, not `test_that()` blocks, so don't try to match an exact number — just confirm `FAIL 0`.)

- [ ] **Step 5: Commit**

```bash
git add R/zone_grouping.R tests/testthat/test-zone-grouping.R
git commit -m "feat: add .validate_zone_assignments for zone grouping input"
```

---

### Task 2: Zone center-finding (medoid) — `.find_zone_centers()`

**Files:**
- Modify: `R/zone_grouping.R` (append)
- Modify: `tests/testthat/test-zone-grouping.R` (append)

**Interfaces:**
- Consumes: `.validate_zone_assignments()` output shape (`tract_id`, `zone_id`, `partition_id`), from Task 1.
- Produces: `.find_zone_centers(zone_assignments, sparse_distances, call = rlang::caller_env())` → tibble with `zone_id`, `partition_id`, `center_tract_id`, `n_tracts`. Aborts, listing every affected zone, when a zone's within-zone tract-pair distances are incomplete.

- [ ] **Step 1: Append failing tests**

```r
test_that(".find_zone_centers picks the medoid minimizing total distance", {
  za <- surveyzones:::.validate_zone_assignments(data.frame(
    tract_id = c("T1", "T2", "T3"),
    zone_id  = c("Z1", "Z1", "Z1")
  ))
  # Positions: T1=0, T2=1, T3=3 -> sums: T1=4, T2=3, T3=5 -> medoid T2
  positions <- c(T1 = 0, T2 = 1, T3 = 3)
  ids <- names(positions)
  dists <- tidyr::expand_grid(origin_id = ids, destination_id = ids) |>
    dplyr::filter(origin_id != destination_id) |>
    dplyr::mutate(distance = abs(positions[origin_id] - positions[destination_id]))

  centers <- surveyzones:::.find_zone_centers(za, dists)
  expect_equal(centers$center_tract_id[centers$zone_id == "Z1"], "T2")
})

test_that(".find_zone_centers handles single-tract zones", {
  za <- surveyzones:::.validate_zone_assignments(data.frame(
    tract_id = "T1", zone_id = "Z1"
  ))
  centers <- surveyzones:::.find_zone_centers(za, tibble::tibble(
    origin_id = character(0), destination_id = character(0), distance = numeric(0)
  ))
  expect_equal(centers$center_tract_id, "T1")
  expect_equal(centers$n_tracts, 1L)
})

test_that(".find_zone_centers uses symmetrized (averaged) distances", {
  # A-B: only d(A,B)=1 given; A-C: d(A,C)=10 and d(C,A)=10;
  # B-C: d(B,C)=1 and d(C,B)=9 -> averaged 5
  # Sums: A = 1+10 = 11; B = 1+5 = 6; C = 10+5 = 15 -> medoid = B
  za <- surveyzones:::.validate_zone_assignments(data.frame(
    tract_id = c("A", "B", "C"), zone_id = c("Z1", "Z1", "Z1")
  ))
  dists <- tibble::tibble(
    origin_id      = c("A", "A", "C", "B", "C"),
    destination_id = c("B", "C", "A", "C", "B"),
    distance       = c(1, 10, 10, 1, 9)
  )
  centers <- surveyzones:::.find_zone_centers(za, dists)
  expect_equal(centers$center_tract_id, "B")
})

test_that(".find_zone_centers breaks ties deterministically (lowest tract_id)", {
  za <- surveyzones:::.validate_zone_assignments(data.frame(
    tract_id = c("T1", "T2", "T3"), zone_id = c("Z1", "Z1", "Z1")
  ))
  dists <- tibble::tibble(
    origin_id      = c("T1", "T1", "T2", "T2", "T3", "T3"),
    destination_id = c("T2", "T3", "T1", "T3", "T1", "T2"),
    distance       = c(2, 2, 2, 2, 2, 2)
  )
  centers <- surveyzones:::.find_zone_centers(za, dists)
  expect_equal(centers$center_tract_id, "T1")
})

test_that(".find_zone_centers aborts on incomplete within-zone distances", {
  za <- surveyzones:::.validate_zone_assignments(data.frame(
    tract_id = c("T1", "T2", "T3"), zone_id = c("Z1", "Z1", "Z1")
  ))
  # T1-T3 pair missing entirely
  dists <- tibble::tibble(
    origin_id      = c("T1", "T2"),
    destination_id = c("T2", "T1"),
    distance       = c(2, 2)
  )
  expect_error(
    surveyzones:::.find_zone_centers(za, dists),
    "incomplete"
  )
})

test_that(".find_zone_centers lists all incomplete zones in one error", {
  za <- surveyzones:::.validate_zone_assignments(data.frame(
    tract_id = c("T1", "T2", "T3", "T4"),
    zone_id  = c("Z1", "Z1", "Z2", "Z2")
  ))
  dists <- tibble::tibble(
    origin_id = character(0), destination_id = character(0), distance = numeric(0)
  )
  err <- tryCatch(
    surveyzones:::.find_zone_centers(za, dists),
    error = function(e) e
  )
  expect_true(grepl("Z1", conditionMessage(err)))
  expect_true(grepl("Z2", conditionMessage(err)))
})
```

- [ ] **Step 2: Run tests to confirm they fail**

```bash
Rscript -e "devtools::load_all(quiet = TRUE); testthat::test_file('tests/testthat/test-zone-grouping.R')"
```

Expected: the 6 new tests fail with `could not find function ".find_zone_centers"`; the 5 Task 1 tests still pass.

- [ ] **Step 3: Append the implementation to `R/zone_grouping.R`**

```r
#' Compute Symmetrized Total Distances Within a Zone
#'
#' @param tract_ids Character vector of tract IDs in one zone (length >= 2).
#' @param sparse_distances Tract-level sparse distance tibble.
#'
#' @return A named numeric vector (sum of symmetrized distances per
#'   candidate tract), or `NULL` if any pair has no distance in either
#'   direction.
#' @keywords internal
.zone_medoid_sums <- function(tract_ids, sparse_distances) {
  n <- length(tract_ids)

  within <- sparse_distances |>
    dplyr::filter(
      .data$origin_id %in% tract_ids,
      .data$destination_id %in% tract_ids,
      .data$origin_id != .data$destination_id
    ) |>
    dplyr::mutate(
      a = pmin(.data$origin_id, .data$destination_id),
      b = pmax(.data$origin_id, .data$destination_id)
    ) |>
    dplyr::summarise(
      distance = mean(.data$distance),
      .by = c("a", "b")
    )

  expected_pairs <- n * (n - 1L) %/% 2L
  if (nrow(within) < expected_pairs) {
    return(NULL)
  }

  sums <- stats::setNames(numeric(n), tract_ids)
  for (i in seq_len(nrow(within))) {
    sums[within$a[i]] <- sums[within$a[i]] + within$distance[i]
    sums[within$b[i]] <- sums[within$b[i]] + within$distance[i]
  }
  sums
}


#' Find Zone Centers by Minimum Total Distance (Medoid)
#'
#' For each zone, finds the tract minimizing the sum of (symmetrized)
#' distances to every other tract in the same zone. Zones with a single
#' tract trivially center on themselves.
#'
#' @param zone_assignments A validated tibble as returned by
#'   [.validate_zone_assignments()] (`tract_id`, `zone_id`, `partition_id`).
#' @param sparse_distances Tract-level sparse distance tibble (`origin_id`,
#'   `destination_id`, `distance`).
#' @param call Caller environment for error reporting.
#'
#' @return A tibble with `zone_id`, `partition_id`, `center_tract_id`,
#'   `n_tracts`.
#' @keywords internal
.find_zone_centers <- function(zone_assignments, sparse_distances,
                               call = rlang::caller_env()) {
  nested <- zone_assignments |>
    tidyr::nest(data = "tract_id") |>
    dplyr::mutate(
      n_tracts = purrr::map_int(.data$data, nrow),
      center_tract_id = purrr::map2_chr(
        .data$data, .data$n_tracts,
        \(d, n) {
          if (n == 1L) return(d$tract_id[1])
          sums <- .zone_medoid_sums(d$tract_id, sparse_distances)
          if (is.null(sums)) return(NA_character_)
          best <- min(sums)
          sort(names(sums)[sums == best])[1]
        }
      )
    )

  incomplete <- nested |>
    dplyr::filter(is.na(.data$center_tract_id)) |>
    dplyr::pull(.data$zone_id)

  if (length(incomplete) > 0) {
    cli::cli_abort(
      c(
        "Cannot find a center for zone{?s} with incomplete tract-pair distances: {.val {incomplete}}.",
        "i" = "Use {.code complete_distances = TRUE} with {.arg access_points} to fill missing pairs."
      ),
      call = call
    )
  }

  nested |>
    dplyr::select("zone_id", "partition_id", "center_tract_id", "n_tracts")
}
```

- [ ] **Step 4: Run tests — expect PASS**

```bash
Rscript -e "devtools::load_all(quiet = TRUE); testthat::test_file('tests/testthat/test-zone-grouping.R')"
```

Expected: `[ FAIL 0 | ... ]` — 0 failures.

- [ ] **Step 5: Commit**

```bash
git add R/zone_grouping.R tests/testthat/test-zone-grouping.R
git commit -m "feat: add .find_zone_centers medoid finder for zone grouping"
```

---

### Task 3: Shared zone-to-zone distance helper — `.zone_dist_sparse()`

**Files:**
- Modify: `R/zone_grouping.R` (append)
- Modify: `tests/testthat/test-zone-grouping.R` (append)

**Interfaces:**
- Consumes: `.build_zone_dist_matrix(zone_ids, assignments, sparse_distances)` and `.matrix_to_sparse_dt(mat, origin_ids, destination_ids, D_max)` — both existing, unmodified, from `R/sequencing.R` and `R/distances.R`.
- Produces: `.zone_dist_sparse(zone_ids, zone_assignments, sparse_distances)` → tibble with `origin_id`, `destination_id`, `distance` (zone IDs as the origin/destination), self-pairs and unreachable (`Inf`) pairs dropped. Used by both `.group_zones_pmedian()` (Task 4) and `.label_zone_groups()` (Task 5).

- [ ] **Step 1: Append failing tests**

```r
test_that(".zone_dist_sparse melts the zone distance matrix without self-pairs", {
  assignments <- tibble::tibble(
    tract_id = paste0("T", 1:4),
    zone_id  = c("Z1", "Z1", "Z2", "Z2"),
    partition_id = "P1"
  )
  positions <- c(T1 = 0, T2 = 1, T3 = 2, T4 = 3)
  ids <- names(positions)
  dists <- tidyr::expand_grid(origin_id = ids, destination_id = ids) |>
    dplyr::filter(origin_id != destination_id) |>
    dplyr::mutate(distance = abs(positions[origin_id] - positions[destination_id]))

  out <- surveyzones:::.zone_dist_sparse(c("Z1", "Z2"), assignments, dists)

  expect_true(all(c("origin_id", "destination_id", "distance") %in% names(out)))
  expect_false(any(out$origin_id == out$destination_id))
  # min pairwise: T2-T3 = 1
  expect_equal(unique(out$distance), 1)
})

test_that(".zone_dist_sparse drops unreachable (Inf) zone pairs", {
  assignments <- tibble::tibble(
    tract_id = c("T1", "T2", "T3"),
    zone_id  = c("Z1", "Z2", "Z3"),
    partition_id = "P1"
  )
  # No distance data involving T3 at all
  dists <- tibble::tibble(
    origin_id = c("T1", "T2"), destination_id = c("T2", "T1"), distance = c(5, 5)
  )
  out <- surveyzones:::.zone_dist_sparse(c("Z1", "Z2", "Z3"), assignments, dists)
  expect_false(any(out$origin_id == "Z3" | out$destination_id == "Z3"))
})
```

- [ ] **Step 2: Run tests to confirm they fail**

```bash
Rscript -e "devtools::load_all(quiet = TRUE); testthat::test_file('tests/testthat/test-zone-grouping.R')"
```

Expected: the 2 new tests fail with `could not find function ".zone_dist_sparse"`; the 11 prior tests still pass.

- [ ] **Step 3: Append the implementation**

```r
#' Build a Sparse Zone-to-Zone Distance Table
#'
#' Wraps [.build_zone_dist_matrix()] (min pairwise tract distance) and
#' melts it into the sparse `origin_id`/`destination_id`/`distance` format
#' used elsewhere in the package, dropping self-pairs and any unreachable
#' (`Inf`) pairs.
#'
#' @param zone_ids Character vector of zone IDs.
#' @param zone_assignments Validated zone assignments (`tract_id`,
#'   `zone_id`).
#' @param sparse_distances Tract-level sparse distance tibble.
#'
#' @return A tibble with `origin_id`, `destination_id`, `distance`
#'   (zone IDs).
#' @keywords internal
.zone_dist_sparse <- function(zone_ids, zone_assignments, sparse_distances) {
  dist_mat <- .build_zone_dist_matrix(zone_ids, zone_assignments, sparse_distances)
  .matrix_to_sparse_dt(
    dist_mat, rownames(dist_mat), colnames(dist_mat), D_max = Inf
  ) |>
    dplyr::filter(is.finite(.data$distance))
}
```

- [ ] **Step 4: Run tests — expect PASS**

```bash
Rscript -e "devtools::load_all(quiet = TRUE); testthat::test_file('tests/testthat/test-zone-grouping.R')"
```

Expected: `[ FAIL 0 | ... ]` — 0 failures.

- [ ] **Step 5: Commit**

```bash
git add R/zone_grouping.R tests/testthat/test-zone-grouping.R
git commit -m "feat: add .zone_dist_sparse shared zone-to-zone distance helper"
```

---

### Task 4: Group zones within D_max — `.group_zones_pmedian()`

**Files:**
- Modify: `R/zone_grouping.R` (append)
- Modify: `tests/testthat/test-zone-grouping.R` (append)

**Interfaces:**
- Consumes: `.zone_dist_sparse()` (Task 3); `surveyzones_build_zones_single()` and `surveyzones_build_zones_single_mem` (existing, from `R/zones.R` / `R/zzz.R`, unmodified).
- Produces: `.group_zones_pmedian(zone_ids, zone_assignments, sparse_distances, D_max, solver, max_time, rel_tol, verbose, use_cache)` → list with `assignments` (tibble: `zone_id`, `group_id`), `groups` (tibble: `group_id`, `center_zone_id`, `n_zones`, `diameter`), `diagnostics` (list: `solver_status`, `objective_value`, `n_variables`, `solve_time`, `optimality_gap`). Input `zone_ids` is scoped to a single partition. Every zone in the output is within `D_max` of its group's `center_zone_id`.

- [ ] **Step 1: Append failing tests**

```r
test_that(".group_zones_pmedian binds close zones into one group", {
  skip_if_not_installed("ROI.plugin.glpk")
  za <- surveyzones:::.validate_zone_assignments(data.frame(
    tract_id = paste0("T", 1:4),
    zone_id  = paste0("Z", 1:4)
  ))
  positions <- c(T1 = 0, T2 = 1, T3 = 2, T4 = 3)
  ids <- names(positions)
  dists <- tidyr::expand_grid(origin_id = ids, destination_id = ids) |>
    dplyr::filter(origin_id != destination_id) |>
    dplyr::mutate(distance = abs(positions[origin_id] - positions[destination_id]))

  out <- surveyzones:::.group_zones_pmedian(
    zone_ids = paste0("Z", 1:4), zone_assignments = za, sparse_distances = dists,
    D_max = 10, solver = "glpk", max_time = 30, rel_tol = 0.01,
    verbose = FALSE, use_cache = FALSE
  )

  expect_equal(nrow(out$groups), 1L)
  expect_equal(out$groups$n_zones, 4L)
})

test_that(".group_zones_pmedian keeps every zone within D_max of its group center", {
  skip_if_not_installed("ROI.plugin.glpk")
  set.seed(1)
  n <- 8
  positions <- stats::runif(n, 0, 50)
  ids <- paste0("Z", seq_len(n))
  names(positions) <- ids
  za <- surveyzones:::.validate_zone_assignments(data.frame(
    tract_id = ids, zone_id = ids
  ))
  dists <- tidyr::expand_grid(origin_id = ids, destination_id = ids) |>
    dplyr::filter(origin_id != destination_id) |>
    dplyr::mutate(distance = abs(positions[origin_id] - positions[destination_id]))

  out <- surveyzones:::.group_zones_pmedian(
    zone_ids = ids, zone_assignments = za, sparse_distances = dists,
    D_max = 10, solver = "glpk", max_time = 30, rel_tol = 0.01,
    verbose = FALSE, use_cache = FALSE
  )

  dist_lookup <- dists |>
    dplyr::mutate(key = paste(origin_id, destination_id)) |>
    dplyr::select(key, distance) |>
    tibble::deframe()

  checked <- out$assignments |>
    dplyr::left_join(out$groups, by = "group_id") |>
    dplyr::filter(.data$zone_id != .data$center_zone_id) |>
    dplyr::mutate(d = dist_lookup[paste(.data$zone_id, .data$center_zone_id)])

  expect_true(all(checked$d <= 10))
})

test_that(".group_zones_pmedian isolates an out-of-range zone as a singleton", {
  skip_if_not_installed("ROI.plugin.glpk")
  ids <- c("Z1", "Z2", "Z3", "Z4")
  positions <- c(Z1 = 0, Z2 = 1, Z3 = 2, Z4 = 1000)
  za <- surveyzones:::.validate_zone_assignments(data.frame(
    tract_id = ids, zone_id = ids
  ))
  dists <- tidyr::expand_grid(origin_id = ids, destination_id = ids) |>
    dplyr::filter(origin_id != destination_id) |>
    dplyr::mutate(distance = abs(positions[origin_id] - positions[destination_id]))

  out <- surveyzones:::.group_zones_pmedian(
    zone_ids = ids, zone_assignments = za, sparse_distances = dists,
    D_max = 10, solver = "glpk", max_time = 30, rel_tol = 0.01,
    verbose = FALSE, use_cache = FALSE
  )

  z4_group <- out$assignments$group_id[out$assignments$zone_id == "Z4"]
  n_zones_in_z4_group <- out$groups$n_zones[out$groups$group_id == z4_group]
  expect_equal(n_zones_in_z4_group, 1L)
})

test_that(".group_zones_pmedian handles a single-zone partition trivially", {
  za <- surveyzones:::.validate_zone_assignments(data.frame(
    tract_id = "T1", zone_id = "Z1"
  ))
  out <- surveyzones:::.group_zones_pmedian(
    zone_ids = "Z1", zone_assignments = za,
    sparse_distances = tibble::tibble(
      origin_id = character(0), destination_id = character(0), distance = numeric(0)
    ),
    D_max = 10, solver = "glpk", max_time = 30, rel_tol = 0.01,
    verbose = FALSE, use_cache = FALSE
  )
  expect_equal(nrow(out$groups), 1L)
  expect_equal(out$groups$group_id, "Z1")
})

test_that(".group_zones_pmedian with solver='spopt' matches glpk's grouping structure", {
  skip_if_not_installed("ROI.plugin.glpk")
  skip_if_not_installed("spopt")
  za <- surveyzones:::.validate_zone_assignments(data.frame(
    tract_id = paste0("T", 1:4),
    zone_id  = paste0("Z", 1:4)
  ))
  positions <- c(T1 = 0, T2 = 1, T3 = 2, T4 = 3)
  ids <- names(positions)
  dists <- tidyr::expand_grid(origin_id = ids, destination_id = ids) |>
    dplyr::filter(origin_id != destination_id) |>
    dplyr::mutate(distance = abs(positions[origin_id] - positions[destination_id]))

  out_glpk <- surveyzones:::.group_zones_pmedian(
    zone_ids = paste0("Z", 1:4), zone_assignments = za, sparse_distances = dists,
    D_max = 10, solver = "glpk", max_time = 30, rel_tol = 0.01,
    verbose = FALSE, use_cache = FALSE
  )
  out_spopt <- surveyzones:::.group_zones_pmedian(
    zone_ids = paste0("Z", 1:4), zone_assignments = za, sparse_distances = dists,
    D_max = 10, solver = "spopt", max_time = 30, rel_tol = 0.01,
    verbose = FALSE, use_cache = FALSE
  )

  expect_equal(nrow(out_glpk$groups), nrow(out_spopt$groups))
  expect_equal(out_glpk$groups$n_zones, out_spopt$groups$n_zones)
})

test_that(".group_zones_pmedian use_cache=TRUE returns identical results on repeat calls", {
  skip_if_not_installed("ROI.plugin.glpk")
  za <- surveyzones:::.validate_zone_assignments(data.frame(
    tract_id = paste0("T", 1:4),
    zone_id  = paste0("Z", 1:4)
  ))
  positions <- c(T1 = 0, T2 = 1, T3 = 2, T4 = 3)
  ids <- names(positions)
  dists <- tidyr::expand_grid(origin_id = ids, destination_id = ids) |>
    dplyr::filter(origin_id != destination_id) |>
    dplyr::mutate(distance = abs(positions[origin_id] - positions[destination_id]))

  args <- list(
    zone_ids = paste0("Z", 1:4), zone_assignments = za, sparse_distances = dists,
    D_max = 10, solver = "glpk", max_time = 30, rel_tol = 0.01,
    verbose = FALSE, use_cache = TRUE
  )
  first <- do.call(surveyzones:::.group_zones_pmedian, args)
  second <- do.call(surveyzones:::.group_zones_pmedian, args)

  expect_equal(first$groups, second$groups)
  expect_equal(first$assignments, second$assignments)
})
```

- [ ] **Step 2: Run tests to confirm they fail**

```bash
Rscript -e "devtools::load_all(quiet = TRUE); testthat::test_file('tests/testthat/test-zone-grouping.R')"
```

Expected: the 6 new tests fail with `could not find function ".group_zones_pmedian"`; the 13 prior tests still pass.

- [ ] **Step 3: Append the implementation**

```r
#' Group Zones Within D_max Using the Reused Tract-Level p-Median Solver
#'
#' Treats zones as the atomic units of a second p-median solve: zone
#' centers become candidate "centers", and the same MILP machinery used
#' for tract-to-zone assignment now bounds every zone to within `D_max`
#' of its group's center zone.
#'
#' @param zone_ids Character vector of zone IDs in a single partition.
#' @param zone_assignments Validated zone assignments (`tract_id`,
#'   `zone_id`, `partition_id`) — used to compute the zone-to-zone
#'   distance matrix via minimum pairwise tract distance.
#' @param sparse_distances Tract-level sparse distance tibble.
#' @param D_max Numeric scalar. Grouping distance bound.
#' @param solver,max_time,rel_tol,verbose,use_cache Passed through to the
#'   reused solver, same meaning as in [surveyzones_build_zones()].
#'
#' @return A list with `assignments` (`zone_id`, `group_id`), `groups`
#'   (`group_id`, `center_zone_id`, `n_zones`, `diameter`), and
#'   `diagnostics`.
#' @keywords internal
.group_zones_pmedian <- function(zone_ids, zone_assignments, sparse_distances,
                                 D_max, solver, max_time, rel_tol, verbose,
                                 use_cache) {
  if (length(zone_ids) == 1L) {
    return(list(
      assignments = tibble::tibble(zone_id = zone_ids, group_id = zone_ids),
      groups = tibble::tibble(
        group_id = zone_ids, center_zone_id = zone_ids,
        n_zones = 1L, diameter = 0
      ),
      diagnostics = list(
        solver_status = "optimal", objective_value = 0,
        n_variables = 0L, solve_time = 0, optimality_gap = 0
      )
    ))
  }

  zone_sparse <- .zone_dist_sparse(zone_ids, zone_assignments, sparse_distances)

  synthetic_tracts <- tibble::tibble(
    tract_id = zone_ids,
    expected_service_time = 1
  )

  solve_single <- if (use_cache) {
    surveyzones_build_zones_single_mem
  } else {
    surveyzones_build_zones_single
  }

  result <- solve_single(
    full_sparse_distances = zone_sparse,
    tracts = synthetic_tracts,
    D_max = D_max,
    max_workload_per_zone = Inf,
    solver = solver,
    max_time = max_time,
    rel_tol = rel_tol,
    verbose = verbose,
    strategy = "auto"
  )

  list(
    assignments = tibble::tibble(
      zone_id = result$assignments$tract_id,
      group_id = result$assignments$zone_id
    ),
    groups = tibble::tibble(
      group_id = result$zones$zone_id,
      center_zone_id = result$zones$center_tract_id,
      n_zones = result$zones$n_tracts,
      diameter = result$zones$diameter
    ),
    diagnostics = result$diagnostics
  )
}
```

- [ ] **Step 4: Run tests — expect PASS**

```bash
Rscript -e "devtools::load_all(quiet = TRUE); testthat::test_file('tests/testthat/test-zone-grouping.R')"
```

Expected: `[ FAIL 0 | ... ]` — 0 failures.

- [ ] **Step 5: Commit**

```bash
git add R/zone_grouping.R tests/testthat/test-zone-grouping.R
git commit -m "feat: add .group_zones_pmedian, reusing the tract-level p-median solver"
```

---

### Task 5: Two independent TSPs and labeling — `.label_zone_groups()`

**Files:**
- Modify: `R/zone_grouping.R` (append)
- Modify: `tests/testthat/test-zone-grouping.R` (append)

**Interfaces:**
- Consumes: `.zone_dist_sparse()` (Task 3); `.seriate_1d(ids, sparse_distances, method = "TSP", control = NULL)` (existing, from `R/sequencing.R`, unmodified — returns `ids` reordered); `sf::st_distance()` (external, for starting-point anchoring, same pattern as `surveyzones_complete_distances()`).
- Produces:
  - `.anchor_partition_order(order_ids, zone_centers, starting_point, access_points)` → `order_ids`, reversed when its *last* element's representative tract is closer to `starting_point` than its *first* element's. A TSP path has exactly two valid global orientations; this picks the one whose first element is closer to the start, without re-solving the TSP with a fixed start node. No-op when `order_ids` has length <= 1 or `starting_point` is `NULL`.
  - `.label_zone_groups(zone_ids, groups, zone_assignments, sparse_distances, zone_centers = NULL, starting_point = NULL, access_points = NULL)` → list with `zone_labels` (tibble: `zone_id`, `zzz`) and `group_labels` (tibble: `group_id`, `ggg`). `groups` is the `groups` tibble from `.group_zones_pmedian()` (Task 4), scoped to one partition. `zzz`/`ggg` are `sprintf("%03d", ...)` strings, each a `1..n` permutation within the partition. The two TSPs are independent — a group's members may be non-contiguous in `zzz`. `zone_centers` is a named character vector (`zone_id` -> representative `tract_id`, scoped to this partition) used only when `starting_point` is supplied, to anchor both TSPs (via `.anchor_partition_order()`) so the first zone and the first group are both the one closest to the partition's starting point.

- [ ] **Step 1: Append failing tests**

```r
test_that(".label_zone_groups produces 3-digit zzz and ggg permutations", {
  ids <- paste0("Z", 1:4)
  assignments <- tibble::tibble(
    tract_id = paste0("T", 1:4), zone_id = ids, partition_id = "P1"
  )
  positions <- c(T1 = 0, T2 = 1, T3 = 2, T4 = 3)
  dists <- tidyr::expand_grid(origin_id = names(positions), destination_id = names(positions)) |>
    dplyr::filter(origin_id != destination_id) |>
    dplyr::mutate(distance = abs(positions[origin_id] - positions[destination_id]))

  groups <- tibble::tibble(
    group_id = c("G1", "G2"),
    center_zone_id = c("Z1", "Z4"),
    n_zones = c(2L, 2L), diameter = c(1, 1)
  )

  out <- surveyzones:::.label_zone_groups(ids, groups, assignments, dists)

  expect_equal(sort(out$zone_labels$zzz), sprintf("%03d", 1:4))
  expect_equal(sort(out$group_labels$ggg), sprintf("%03d", 1:2))
})

test_that(".label_zone_groups keeps zzz independent of group membership", {
  # 4 zones on a line: Z1(0) Z2(1) Z3(2) Z4(3)
  # Groups: G1 = {Z1, Z3}, G2 = {Z2, Z4} -- interleaved on the line
  ids <- paste0("Z", 1:4)
  assignments <- tibble::tibble(tract_id = ids, zone_id = ids, partition_id = "P1")
  positions <- c(Z1 = 0, Z2 = 1, Z3 = 2, Z4 = 3)
  dists <- tidyr::expand_grid(origin_id = ids, destination_id = ids) |>
    dplyr::filter(origin_id != destination_id) |>
    dplyr::mutate(distance = abs(positions[origin_id] - positions[destination_id]))

  groups <- tibble::tibble(
    group_id = c("G1", "G2"), center_zone_id = c("Z1", "Z2"),
    n_zones = c(2L, 2L), diameter = c(2, 2)
  )

  out <- surveyzones:::.label_zone_groups(ids, groups, assignments, dists)

  zzz_num <- stats::setNames(as.integer(out$zone_labels$zzz), out$zone_labels$zone_id)

  # The zone TSP recovers the spatial line order (Z1..Z4 or its reverse)
  ordered_by_zzz <- names(sort(zzz_num))
  in_position_order <- all(ordered_by_zzz == paste0("Z", 1:4)) ||
    all(ordered_by_zzz == paste0("Z", 4:1))
  expect_true(in_position_order)

  # G1 = {Z1, Z3}: their zzz positions are 2 apart, not adjacent -- proving
  # the zone TSP was not constrained by group membership
  expect_equal(abs(zzz_num[["Z1"]] - zzz_num[["Z3"]]), 2L)
})

test_that(".label_zone_groups handles a single zone and single group", {
  ids <- "Z1"
  assignments <- tibble::tibble(tract_id = "T1", zone_id = "Z1", partition_id = "P1")
  dists <- tibble::tibble(
    origin_id = character(0), destination_id = character(0), distance = numeric(0)
  )
  groups <- tibble::tibble(
    group_id = "Z1", center_zone_id = "Z1", n_zones = 1L, diameter = 0
  )

  out <- surveyzones:::.label_zone_groups(ids, groups, assignments, dists)
  expect_equal(out$zone_labels$zzz, "001")
  expect_equal(out$group_labels$ggg, "001")
})

test_that(".anchor_partition_order reverses order when the last element is closer to start", {
  ids <- c("Z1", "Z2", "Z3")
  zone_centers <- c(Z1 = "T1", Z2 = "T2", Z3 = "T3")
  access_points <- sf::st_as_sf(
    data.frame(
      tract_id = c("T1", "T2", "T3"),
      lon = c(0, 1, 2), lat = c(0, 0, 0)
    ),
    coords = c("lon", "lat"), crs = 4326
  )
  starting_point <- sf::st_as_sf(
    data.frame(lon = 3, lat = 0), coords = c("lon", "lat"), crs = 4326
  )

  out <- surveyzones:::.anchor_partition_order(ids, zone_centers, starting_point, access_points)
  expect_equal(out, c("Z3", "Z2", "Z1"))
})

test_that(".anchor_partition_order keeps order when the first element is already closer", {
  ids <- c("Z1", "Z2", "Z3")
  zone_centers <- c(Z1 = "T1", Z2 = "T2", Z3 = "T3")
  access_points <- sf::st_as_sf(
    data.frame(
      tract_id = c("T1", "T2", "T3"),
      lon = c(0, 1, 2), lat = c(0, 0, 0)
    ),
    coords = c("lon", "lat"), crs = 4326
  )
  starting_point <- sf::st_as_sf(
    data.frame(lon = -1, lat = 0), coords = c("lon", "lat"), crs = 4326
  )

  out <- surveyzones:::.anchor_partition_order(ids, zone_centers, starting_point, access_points)
  expect_equal(out, c("Z1", "Z2", "Z3"))
})

test_that(".anchor_partition_order is a no-op for a single element or NULL starting_point", {
  expect_equal(
    surveyzones:::.anchor_partition_order("Z1", c(Z1 = "T1"), NULL, NULL),
    "Z1"
  )
  expect_equal(
    surveyzones:::.anchor_partition_order(c("Z1", "Z2"), c(Z1 = "T1", Z2 = "T2"), NULL, NULL),
    c("Z1", "Z2")
  )
})

test_that(".label_zone_groups anchors both zzz and ggg to the partition's starting point", {
  # 4 zones on a line: Z1(0) Z2(1) Z3(2) Z4(3); groups G1={Z1,Z2}, G2={Z3,Z4}.
  # Starting point is right at Z4 -> zzz should start at Z4 and ggg should
  # start with G2 (whose center is Z4).
  ids <- paste0("Z", 1:4)
  assignments <- tibble::tibble(tract_id = ids, zone_id = ids, partition_id = "P1")
  positions <- c(Z1 = 0, Z2 = 1, Z3 = 2, Z4 = 3)
  dists <- tidyr::expand_grid(origin_id = ids, destination_id = ids) |>
    dplyr::filter(origin_id != destination_id) |>
    dplyr::mutate(distance = abs(positions[origin_id] - positions[destination_id]))

  groups <- tibble::tibble(
    group_id = c("G1", "G2"), center_zone_id = c("Z1", "Z4"),
    n_zones = c(2L, 2L), diameter = c(1, 1)
  )

  zone_centers <- c(Z1 = "Z1", Z2 = "Z2", Z3 = "Z3", Z4 = "Z4")
  access_points <- sf::st_as_sf(
    data.frame(tract_id = ids, lon = unname(positions), lat = 0),
    coords = c("lon", "lat"), crs = 4326
  )
  starting_point <- sf::st_as_sf(
    data.frame(lon = 3, lat = 0), coords = c("lon", "lat"), crs = 4326
  )

  out <- surveyzones:::.label_zone_groups(
    ids, groups, assignments, dists,
    zone_centers = zone_centers, starting_point = starting_point,
    access_points = access_points
  )

  expect_equal(out$zone_labels$zone_id[out$zone_labels$zzz == "001"], "Z4")
  expect_equal(out$group_labels$group_id[out$group_labels$ggg == "001"], "G2")
})
```

- [ ] **Step 2: Run tests to confirm they fail**

```bash
Rscript -e "devtools::load_all(quiet = TRUE); testthat::test_file('tests/testthat/test-zone-grouping.R')"
```

Expected: the 7 new tests fail with `could not find function ".label_zone_groups"` (and `.anchor_partition_order`); all prior tests still pass.

- [ ] **Step 3: Append the implementation**

```r
#' Reorient a Zone/Group Visit Order Toward a Partition's Starting Point
#'
#' Compares the endpoint tracts of an ordered ID sequence (each ID mapped
#' to its representative tract) to a partition's starting-point
#' coordinate, and reverses the sequence when its *last* element is
#' closer to the starting point than its *first* element -- so the
#' sequence starts at the (endpoint-optimal) closer end. TSP paths have
#' exactly two valid global orientations; this picks the better one
#' without re-solving the TSP with a fixed start node.
#'
#' @param order_ids Character vector, already TSP-ordered (zone IDs or
#'   group center-zone IDs).
#' @param zone_centers Named character vector: representative `tract_id`
#'   for each ID in `order_ids` (names = the IDs).
#' @param starting_point A one-row sf POINT object (the partition's
#'   starting point), or `NULL` to skip anchoring.
#' @param access_points An sf object with POINT geometries and a
#'   `tract_id` column, covering every tract in `zone_centers`.
#'
#' @return `order_ids`, possibly reversed.
#' @keywords internal
.anchor_partition_order <- function(order_ids, zone_centers, starting_point,
                                    access_points) {
  if (length(order_ids) <= 1L || is.null(starting_point)) return(order_ids)

  first_tract <- zone_centers[[order_ids[1]]]
  last_tract <- zone_centers[[order_ids[length(order_ids)]]]

  first_pt <- access_points[access_points$tract_id == first_tract, ]
  last_pt <- access_points[access_points$tract_id == last_tract, ]

  d_first <- as.numeric(sf::st_distance(starting_point, first_pt))
  d_last <- as.numeric(sf::st_distance(starting_point, last_pt))

  if (d_last < d_first) rev(order_ids) else order_ids
}


#' Assign ggg.zzz Labels via Two Independent Partition-Wide TSPs
#'
#' Seriates group centers (for `ggg`) and all zones (for `zzz`)
#' independently within one partition -- the two orders do not need to
#' agree, so a group's members may be non-contiguous in the `zzz`
#' sequence. When a starting point is supplied, both orders are anchored
#' via [.anchor_partition_order()] so the first zone and first group are
#' each the one closest to it.
#'
#' @param zone_ids Character vector of all zone IDs in one partition.
#' @param groups Tibble for this partition with `group_id`,
#'   `center_zone_id` (as returned by [.group_zones_pmedian()]).
#' @param zone_assignments Validated zone assignments.
#' @param sparse_distances Tract-level sparse distance tibble.
#' @param zone_centers Named character vector (`zone_id` -> representative
#'   `tract_id`, scoped to this partition). Only used when
#'   `starting_point` is supplied.
#' @param starting_point A one-row sf POINT object (this partition's
#'   starting point), or `NULL` (default) to skip anchoring.
#' @param access_points An sf object with POINT geometries and a
#'   `tract_id` column. Required when `starting_point` is supplied.
#'
#' @return A list with `zone_labels` (`zone_id`, `zzz`) and
#'   `group_labels` (`group_id`, `ggg`).
#' @keywords internal
.label_zone_groups <- function(zone_ids, groups, zone_assignments,
                               sparse_distances, zone_centers = NULL,
                               starting_point = NULL, access_points = NULL) {
  zone_dists <- .zone_dist_sparse(zone_ids, zone_assignments, sparse_distances)

  zone_order <- .seriate_1d(zone_ids, zone_dists)
  group_order <- .seriate_1d(groups$center_zone_id, zone_dists)

  zone_order <- .anchor_partition_order(
    zone_order, zone_centers, starting_point, access_points
  )
  group_order <- .anchor_partition_order(
    group_order, zone_centers, starting_point, access_points
  )

  zone_labels <- tibble::tibble(
    zone_id = zone_order,
    zzz = sprintf("%03d", seq_along(zone_order))
  )

  group_labels <- tibble::tibble(
    center_zone_id = group_order,
    ggg = sprintf("%03d", seq_along(group_order))
  ) |>
    dplyr::left_join(
      groups |> dplyr::select("group_id", "center_zone_id"),
      by = "center_zone_id"
    ) |>
    dplyr::select("group_id", "ggg")

  list(zone_labels = zone_labels, group_labels = group_labels)
}
```

- [ ] **Step 4: Run tests — expect PASS**

```bash
Rscript -e "devtools::load_all(quiet = TRUE); testthat::test_file('tests/testthat/test-zone-grouping.R')"
```

Expected: `[ FAIL 0 | ... ]` — 0 failures.

- [ ] **Step 5: Commit**

```bash
git add R/zone_grouping.R tests/testthat/test-zone-grouping.R
git commit -m "feat: add .label_zone_groups with starting-point anchoring for ggg/zzz"
```

---

### Task 6: Tract-level sequencing and cross-zone orientation — `.sequence_group_tracts()`

**Files:**
- Modify: `R/zone_grouping.R` (append)
- Modify: `tests/testthat/test-zone-grouping.R` (append)

**Interfaces:**
- Consumes: `surveyzones_sequence_tracts(tract_ids, sparse_distances, method = "TSP", control = NULL)` (existing, exported, from `R/sequencing.R`, unmodified).
- Produces: `.sequence_group_tracts(zone_order, zone_assignments, sparse_distances, method = "TSP", control = NULL)` → tibble with `zone_id`, `tract_id`, `visit_order`. `zone_order` is a character vector of zone IDs already in walk order (the `zzz` TSP order from Task 5, for one partition). Each zone's tract list is oriented so its entry tract is closest to the previous zone's exit tract.

- [ ] **Step 1: Append failing tests**

```r
test_that(".sequence_group_tracts produces a permutation of each zone's tracts", {
  za <- tibble::tibble(
    tract_id = paste0("T", 1:4),
    zone_id  = c("ZA", "ZA", "ZB", "ZB")
  )
  positions <- c(T1 = 0, T2 = 1, T3 = 10, T4 = 11)
  ids <- names(positions)
  dists <- tidyr::expand_grid(origin_id = ids, destination_id = ids) |>
    dplyr::filter(origin_id != destination_id) |>
    dplyr::mutate(distance = abs(positions[origin_id] - positions[destination_id]))

  out <- surveyzones:::.sequence_group_tracts(c("ZA", "ZB"), za, dists)

  for (zid in c("ZA", "ZB")) {
    assigned <- sort(za$tract_id[za$zone_id == zid])
    sequenced <- sort(out$tract_id[out$zone_id == zid])
    expect_equal(assigned, sequenced)
  }
})

test_that(".sequence_group_tracts orients the second zone toward the first zone's exit", {
  za <- tibble::tibble(
    tract_id = paste0("T", 1:6),
    zone_id  = c("ZA", "ZA", "ZA", "ZB", "ZB", "ZB")
  )
  positions <- c(T1 = 0, T2 = 1, T3 = 2, T4 = 10, T5 = 11, T6 = 12)
  ids <- names(positions)
  dists <- tidyr::expand_grid(origin_id = ids, destination_id = ids) |>
    dplyr::filter(origin_id != destination_id) |>
    dplyr::mutate(distance = abs(positions[origin_id] - positions[destination_id]))

  out <- surveyzones:::.sequence_group_tracts(c("ZA", "ZB"), za, dists)

  za_tracts <- out |> dplyr::filter(zone_id == "ZA") |> dplyr::arrange(visit_order)
  zb_tracts <- out |> dplyr::filter(zone_id == "ZB") |> dplyr::arrange(visit_order)

  exit_pos <- positions[za_tracts$tract_id[nrow(za_tracts)]]
  entry_pos <- positions[zb_tracts$tract_id[1]]
  tail_pos <- positions[zb_tracts$tract_id[nrow(zb_tracts)]]

  expect_true(abs(entry_pos - exit_pos) <= abs(tail_pos - exit_pos))
})

test_that(".sequence_group_tracts handles single-tract zones", {
  za <- tibble::tibble(tract_id = c("T1", "T2"), zone_id = c("Z1", "Z2"))
  dists <- tibble::tibble(
    origin_id = character(0), destination_id = character(0), distance = numeric(0)
  )
  out <- surveyzones:::.sequence_group_tracts(c("Z1", "Z2"), za, dists)
  expect_equal(nrow(out), 2L)
  expect_equal(out$visit_order, c(1L, 1L))
})
```

- [ ] **Step 2: Run tests to confirm they fail**

```bash
Rscript -e "devtools::load_all(quiet = TRUE); testthat::test_file('tests/testthat/test-zone-grouping.R')"
```

Expected: the 3 new tests fail with `could not find function ".sequence_group_tracts"`; the 26 prior tests still pass.

- [ ] **Step 3: Append the implementation**

```r
#' Sequence and Orient Tracts Along a Given Zone Visit Order
#'
#' Walks zones in the supplied order, seriates each zone's tracts, and
#' orients each zone's tract sequence so its entry tract is closest to the
#' previous zone's exit tract. Mirrors
#' `.sequence_tracts_with_orientation()` in `R/sequencing.R`, but driven
#' by an externally supplied zone order (here, the `zzz` TSP order from
#' [.label_zone_groups()]) instead of a group+within-group order.
#'
#' @param zone_order Character vector of zone IDs, in visit order
#'   (one partition).
#' @param zone_assignments Validated zone assignments (`tract_id`,
#'   `zone_id`).
#' @param sparse_distances Tract-level sparse distance tibble.
#' @param method,control Passed to [surveyzones_sequence_tracts()].
#'
#' @return A tibble with `zone_id`, `tract_id`, `visit_order`.
#' @keywords internal
.sequence_group_tracts <- function(zone_order, zone_assignments,
                                   sparse_distances, method = "TSP",
                                   control = NULL) {
  dist_key <- paste(
    sparse_distances$origin_id, sparse_distances$destination_id, sep = "|"
  )
  dist_val <- stats::setNames(sparse_distances$distance, dist_key)

  tracts_by_zone <- split(
    zone_assignments$tract_id, zone_assignments$zone_id
  )

  sequences <- vector("list", length(zone_order))
  prev_exit <- NA_character_

  for (j in seq_along(zone_order)) {
    zid <- zone_order[j]
    tract_ids <- tracts_by_zone[[zid]]

    ordered_tracts <- surveyzones_sequence_tracts(
      tract_ids, sparse_distances, method = method, control = control
    )

    if (!is.na(prev_exit) && length(ordered_tracts) > 1) {
      first <- ordered_tracts[1]
      last <- ordered_tracts[length(ordered_tracts)]
      d_first <- dist_val[paste(prev_exit, first, sep = "|")]
      d_last <- dist_val[paste(prev_exit, last, sep = "|")]
      if (!is.na(d_first) && !is.na(d_last) && d_last < d_first) {
        ordered_tracts <- rev(ordered_tracts)
      }
    }

    prev_exit <- ordered_tracts[length(ordered_tracts)]

    sequences[[j]] <- tibble::tibble(
      zone_id = zid,
      tract_id = ordered_tracts,
      visit_order = seq_along(ordered_tracts)
    )
  }

  dplyr::bind_rows(sequences)
}
```

- [ ] **Step 4: Run tests — expect PASS**

```bash
Rscript -e "devtools::load_all(quiet = TRUE); testthat::test_file('tests/testthat/test-zone-grouping.R')"
```

Expected: `[ FAIL 0 | ... ]` — 0 failures.

- [ ] **Step 5: Commit**

```bash
git add R/zone_grouping.R tests/testthat/test-zone-grouping.R
git commit -m "feat: add .sequence_group_tracts for cross-zone tract orientation"
```

---

### Task 7: `surveyzones_zone_groups` S3 class — constructor and print method

**Files:**
- Create: `R/zone_groups_plan.R`
- Modify: `tests/testthat/test-zone-grouping.R` (append)

**Interfaces:**
- Produces: `new_surveyzones_zone_groups(zones, groups, sequence, parameters, diagnostics)` → object of class `surveyzones_zone_groups` (a list with those five named elements). `print.surveyzones_zone_groups()` S3 method.

- [ ] **Step 1: Append failing test**

```r
test_that("new_surveyzones_zone_groups and print work", {
  obj <- surveyzones:::new_surveyzones_zone_groups(
    zones = tibble::tibble(
      zone_id = "Z1", partition_id = "P1", center_tract_id = "T1",
      n_tracts = 1L, group_id = "Z1", ggg = "001", zzz = "001", label = "001.001"
    ),
    groups = tibble::tibble(
      group_id = "Z1", partition_id = "P1", center_zone_id = "Z1",
      n_zones = 1L, diameter = 0, ggg = "001"
    ),
    sequence = tibble::tibble(zone_id = "Z1", tract_id = "T1", visit_order = 1L),
    parameters = list(D_max = 10, solver = "glpk", timestamp = Sys.time()),
    diagnostics = list(solver_status = "optimal")
  )
  expect_s3_class(obj, "surveyzones_zone_groups")
  expect_no_error(print(obj))
})
```

- [ ] **Step 2: Run test to confirm it fails**

```bash
Rscript -e "devtools::load_all(quiet = TRUE); testthat::test_file('tests/testthat/test-zone-grouping.R')"
```

Expected: fails with `could not find function "new_surveyzones_zone_groups"`; the 29 prior tests still pass.

- [ ] **Step 3: Create `R/zone_groups_plan.R`**

```r
#' Create a surveyzones_zone_groups Object
#'
#' Constructor for the S3 class returned by [surveyzones_group_zones()].
#'
#' @param zones A tibble with columns `zone_id`, `partition_id`,
#'   `center_tract_id`, `n_tracts`, `group_id`, `ggg`, `zzz`, `label`.
#' @param groups A tibble with columns `group_id`, `partition_id`,
#'   `center_zone_id`, `n_zones`, `diameter`, `ggg`.
#' @param sequence A tibble with columns `zone_id`, `tract_id`,
#'   `visit_order`.
#' @param parameters A named list of parameters (`D_max`, `solver`,
#'   `timestamp`).
#' @param diagnostics A named list of solver diagnostics.
#'
#' @return A `surveyzones_zone_groups` object.
#' @keywords internal
new_surveyzones_zone_groups <- function(zones, groups, sequence,
                                        parameters, diagnostics) {
  structure(
    list(
      zones = zones,
      groups = groups,
      sequence = sequence,
      parameters = parameters,
      diagnostics = diagnostics
    ),
    class = "surveyzones_zone_groups"
  )
}


#' @export
print.surveyzones_zone_groups <- function(x, ...) {
  n_partitions <- length(unique(x$zones$partition_id))
  n_zones <- nrow(x$zones)
  n_groups <- nrow(x$groups)

  cli::cli_h1("surveyzones zone groups")
  cli::cli_bullets(c(
    "*" = "{n_partitions} partition{?s}, {n_zones} zone{?s}, {n_groups} group{?s}",
    "*" = "D_max = {x$parameters$D_max}"
  ))

  if (n_groups > 0) {
    dm <- x$groups$diameter
    cli::cli_bullets(c(
      "*" = "Group diameter: {round(min(dm, na.rm = TRUE), 2)}–{round(max(dm, na.rm = TRUE), 2)} (mean {round(mean(dm, na.rm = TRUE), 2)})"
    ))
  }

  invisible(x)
}
```

- [ ] **Step 4: Run test — expect PASS**

```bash
Rscript -e "devtools::load_all(quiet = TRUE); testthat::test_file('tests/testthat/test-zone-grouping.R')"
```

Expected: `[ FAIL 0 | ... ]` — 0 failures.

- [ ] **Step 5: Commit**

```bash
git add R/zone_groups_plan.R tests/testthat/test-zone-grouping.R
git commit -m "feat: add surveyzones_zone_groups S3 constructor and print method"
```

---

### Task 8: Top-level orchestration — `surveyzones_group_zones()`

**Files:**
- Modify: `R/zone_grouping.R` (append)
- Modify: `tests/testthat/test-zone-grouping.R` (append)

**Interfaces:**
- Consumes: every function from Tasks 1-7 (`.validate_zone_assignments`, `.find_zone_centers`, `.group_zones_pmedian`, `.label_zone_groups`, `.anchor_partition_order`, `.sequence_group_tracts`, `new_surveyzones_zone_groups`), plus existing `validate_solver()` and `surveyzones_complete_distances()` (unmodified).
- Produces: exported `surveyzones_group_zones(zone_assignments, sparse_distances, D_max, solver = "glpk", max_time = 300, rel_tol = 0.01, verbose = FALSE, use_cache = TRUE, access_points = NULL, speed_kmh = 0.1, complete_distances = TRUE, starting_points = NULL, method = "TSP", control = NULL)` → `surveyzones_zone_groups` object. `starting_points` is an optional sf POINT object with a `partition_id` column, one row per partition; when a partition has a matching row, that partition's `zzz` and `ggg` orders are anchored (via `.label_zone_groups()`'s `starting_point`/`zone_centers` args) so the first zone and first group are the ones closest to it. Partitions absent from `starting_points` keep the unanchored TSP orientation. Requires `access_points` whenever `starting_points` is supplied, regardless of `complete_distances`.

- [ ] **Step 1: Append failing tests**

```r
test_that("surveyzones_group_zones returns a surveyzones_zone_groups object", {
  skip_if_not_installed("ROI.plugin.glpk")
  za <- data.frame(
    tract_id = paste0("T", 1:4),
    zone_id  = c("Z1", "Z1", "Z2", "Z2")
  )
  positions <- c(T1 = 0, T2 = 1, T3 = 2, T4 = 3)
  ids <- names(positions)
  dists <- tidyr::expand_grid(origin_id = ids, destination_id = ids) |>
    dplyr::filter(origin_id != destination_id) |>
    dplyr::mutate(distance = abs(positions[origin_id] - positions[destination_id]))

  result <- surveyzones_group_zones(
    zone_assignments = za,
    sparse_distances = dists,
    D_max = 10,
    solver = "glpk",
    use_cache = FALSE,
    complete_distances = FALSE
  )

  expect_s3_class(result, "surveyzones_zone_groups")
  expect_equal(nrow(result$zones), 2L)
  expect_equal(sort(result$sequence$tract_id), sort(za$tract_id))
})

test_that("surveyzones_group_zones requires access_points when complete_distances = TRUE", {
  za <- data.frame(tract_id = "T1", zone_id = "Z1")
  dists <- tibble::tibble(
    origin_id = character(0), destination_id = character(0), distance = numeric(0)
  )
  expect_error(
    surveyzones_group_zones(za, dists, D_max = 10),
    "access_points"
  )
})

test_that("surveyzones_group_zones resolves centers after haversine completion", {
  skip_if_not_installed("ROI.plugin.glpk")
  za <- data.frame(
    tract_id = c("T1", "T2", "T3"),
    zone_id  = c("Z1", "Z1", "Z1")
  )
  pts <- sf::st_as_sf(
    data.frame(
      tract_id = c("T1", "T2", "T3"),
      lon = c(-38.50, -38.501, -38.502),
      lat = c(-13.00, -13.001, -13.002)
    ),
    coords = c("lon", "lat"), crs = 4326
  )
  # Only one pair supplied; the rest must come from haversine completion
  dists <- tibble::tibble(
    origin_id = c("T1", "T2"), destination_id = c("T2", "T1"), distance = c(1, 1)
  )

  result <- surveyzones_group_zones(
    zone_assignments = za,
    sparse_distances = dists,
    D_max = 100,
    solver = "glpk",
    use_cache = FALSE,
    access_points = pts,
    complete_distances = TRUE
  )

  expect_equal(nrow(result$zones), 1L)
  expect_true(result$zones$center_tract_id %in% za$tract_id)
})

test_that("surveyzones_group_zones anchors the first zone/group per partition to starting_points", {
  skip_if_not_installed("ROI.plugin.glpk")
  # 4 zones on a line per partition: Z1(0) Z2(1) Z3(2) Z4(3); D_max=10 -> 1 group.
  # Starting point at position 3 -> the zone/group nearest it (Z4) must be first.
  za <- tibble::tibble(
    zone_id = paste0("Z", 1:4),
    tract_id = paste0("Z", 1:4),
    partition_id = "P1"
  )
  positions <- c(Z1 = 0, Z2 = 1, Z3 = 2, Z4 = 3)
  pts <- sf::st_as_sf(
    data.frame(tract_id = names(positions), lon = unname(positions), lat = 0),
    coords = c("lon", "lat"), crs = 4326
  )
  dists <- tidyr::expand_grid(origin_id = names(positions), destination_id = names(positions)) |>
    dplyr::filter(origin_id != destination_id) |>
    dplyr::mutate(distance = abs(positions[origin_id] - positions[destination_id]))

  starting_points <- sf::st_as_sf(
    data.frame(partition_id = "P1", lon = 3, lat = 0),
    coords = c("lon", "lat"), crs = 4326
  )

  result <- surveyzones_group_zones(
    zone_assignments = za,
    sparse_distances = dists,
    D_max = 10,
    solver = "glpk",
    use_cache = FALSE,
    access_points = pts,
    complete_distances = FALSE,
    starting_points = starting_points
  )

  expect_equal(result$zones$zone_id[result$zones$zzz == "001"], "Z4")
  expect_equal(result$groups$center_zone_id[result$groups$ggg == "001"], "Z4")
})
```

- [ ] **Step 2: Run tests to confirm they fail**

```bash
Rscript -e "devtools::load_all(quiet = TRUE); testthat::test_file('tests/testthat/test-zone-grouping.R')"
```

Expected: the 4 new tests fail with `could not find function "surveyzones_group_zones"`; the 30 prior tests still pass.

- [ ] **Step 3: Append the implementation**

```r
#' Group Externally-Defined Zones Within a Distance Bound and Label Them
#'
#' Given zones that already exist (not built by [surveyzones_build_zones()]),
#' finds a center for each zone, groups zones within `D_max` of a group
#' center by reusing the same p-median machinery one level up, and assigns
#' each zone a stable `ggg.zzz` label from two independent partition-wide
#' TSPs. Also produces a tract-level visit sequence oriented across zone
#' boundaries.
#'
#' @param zone_assignments A data.frame with `tract_id`, `zone_id`, and
#'   optionally `partition_id`. Every tract in a given `zone_id` must
#'   share the same `partition_id`.
#' @param sparse_distances Tract-level sparse distance tibble (`origin_id`,
#'   `destination_id`, `distance`).
#' @param D_max Numeric scalar. Every zone in a group is within `D_max` of
#'   that group's center zone.
#' @param solver Character scalar. One of `"glpk"` (default), `"highs"`,
#'   `"cbc"`, `"symphony"`, `"spopt"` -- reused for the zone-level p-median.
#' @param max_time,rel_tol,verbose,use_cache Passed through to the reused
#'   solver, same meaning as in [surveyzones_build_zones()].
#' @param access_points An sf object with POINT geometries and a
#'   `tract_id` column. Used to fill missing pairwise distances via
#'   haversine. Required when `complete_distances = TRUE`.
#' @param speed_kmh Numeric scalar. Assumed travel speed for haversine
#'   fill-in (km/h). Default `0.1`.
#' @param complete_distances Logical scalar. When `TRUE` (default),
#'   missing pairwise distances are filled using haversine via
#'   [surveyzones_complete_distances()].
#' @param starting_points An optional sf POINT object with a
#'   `partition_id` column, one row per partition. When a partition has a
#'   matching row, that partition's `zzz` and `ggg` orders are anchored so
#'   the first zone and first group are the ones closest to it (see
#'   [.anchor_partition_order()]). Partitions absent from `starting_points`
#'   keep the unanchored TSP orientation. Requires `access_points`
#'   regardless of `complete_distances`.
#' @param method,control Seriation method/control passed to
#'   [surveyzones_sequence_tracts()] for tract-level ordering. Zone- and
#'   group-level TSPs always use `"TSP"`.
#'
#' @return A `surveyzones_zone_groups` object.
#'
#' @export
surveyzones_group_zones <- function(
  zone_assignments,
  sparse_distances,
  D_max,
  solver = "glpk",
  max_time = 300,
  rel_tol = 0.01,
  verbose = FALSE,
  use_cache = TRUE,
  access_points = NULL,
  speed_kmh = 0.1,
  complete_distances = TRUE,
  starting_points = NULL,
  method = "TSP",
  control = NULL
) {
  D_max <- as.numeric(D_max)
  validate_solver(solver)

  za <- .validate_zone_assignments(zone_assignments)

  if (!is.null(starting_points) && is.null(access_points)) {
    cli::cli_abort(c(
      "{.arg access_points} is required when {.arg starting_points} is supplied.",
      "i" = "Starting-point anchoring needs real coordinates to compute distances."
    ))
  }
  if (!is.null(starting_points) && !"partition_id" %in% names(starting_points)) {
    cli::cli_abort("{.arg starting_points} must contain a {.field partition_id} column.")
  }

  if (complete_distances) {
    if (is.null(access_points)) {
      cli::cli_abort(c(
        "{.arg access_points} is required when {.arg complete_distances} is {.val TRUE}.",
        "i" = "Either pass it explicitly or set {.code complete_distances = FALSE}."
      ))
    }
    sparse_distances <- surveyzones_complete_distances(
      sparse_distances, access_points, speed_kmh = speed_kmh
    )
  }

  centers <- .find_zone_centers(za, sparse_distances)

  partitions <- split(centers$zone_id, centers$partition_id)

  zones_out <- vector("list", length(partitions))
  groups_out <- vector("list", length(partitions))
  sequence_out <- vector("list", length(partitions))
  diagnostics_out <- vector("list", length(partitions))

  for (idx in seq_along(partitions)) {
    pid <- names(partitions)[idx]
    zone_ids <- partitions[[idx]]
    za_part <- za |> dplyr::filter(.data$partition_id == pid)
    centers_part <- centers |> dplyr::filter(.data$partition_id == pid)
    zone_centers_part <- stats::setNames(
      centers_part$center_tract_id, centers_part$zone_id
    )

    starting_point_part <- if (is.null(starting_points)) {
      NULL
    } else {
      sp <- starting_points[starting_points$partition_id == pid, ]
      if (nrow(sp) == 0) NULL else sp[1, ]
    }

    grouped <- .group_zones_pmedian(
      zone_ids = zone_ids, zone_assignments = za_part,
      sparse_distances = sparse_distances, D_max = D_max,
      solver = solver, max_time = max_time, rel_tol = rel_tol,
      verbose = verbose, use_cache = use_cache
    )

    labels <- .label_zone_groups(
      zone_ids = zone_ids, groups = grouped$groups,
      zone_assignments = za_part, sparse_distances = sparse_distances,
      zone_centers = zone_centers_part, starting_point = starting_point_part,
      access_points = access_points
    )

    zone_order <- labels$zone_labels |>
      dplyr::arrange(as.integer(.data$zzz)) |>
      dplyr::pull(.data$zone_id)

    sequence_out[[idx]] <- .sequence_group_tracts(
      zone_order, za_part, sparse_distances, method = method, control = control
    )

    part_zones <- centers |>
      dplyr::filter(.data$partition_id == pid) |>
      dplyr::left_join(grouped$assignments, by = "zone_id") |>
      dplyr::left_join(labels$zone_labels, by = "zone_id") |>
      dplyr::left_join(labels$group_labels, by = "group_id") |>
      dplyr::mutate(label = paste0(.data$ggg, ".", .data$zzz))

    zones_out[[idx]] <- part_zones |>
      dplyr::select(
        "zone_id", "partition_id", "center_tract_id", "n_tracts",
        "group_id", "ggg", "zzz", "label"
      )

    groups_out[[idx]] <- grouped$groups |>
      dplyr::left_join(labels$group_labels, by = "group_id") |>
      dplyr::mutate(partition_id = pid) |>
      dplyr::select(
        "group_id", "partition_id", "center_zone_id", "n_zones",
        "diameter", "ggg"
      )

    diagnostics_out[[idx]] <- grouped$diagnostics
  }

  new_surveyzones_zone_groups(
    zones = dplyr::bind_rows(zones_out),
    groups = dplyr::bind_rows(groups_out),
    sequence = dplyr::bind_rows(sequence_out),
    parameters = list(D_max = D_max, solver = solver, timestamp = Sys.time()),
    diagnostics = list(
      solver_status = purrr::map_chr(diagnostics_out, "solver_status"),
      objective_value = purrr::map_dbl(diagnostics_out, "objective_value"),
      n_variables = purrr::map_int(diagnostics_out, "n_variables"),
      solve_time = purrr::map_dbl(diagnostics_out, "solve_time"),
      optimality_gap = purrr::map_dbl(diagnostics_out, "optimality_gap")
    )
  )
}
```

- [ ] **Step 4: Run tests — expect PASS**

```bash
Rscript -e "devtools::load_all(quiet = TRUE); testthat::test_file('tests/testthat/test-zone-grouping.R')"
```

Expected: `[ FAIL 0 | ... ]` — 0 failures.

- [ ] **Step 5: Commit**

```bash
git add R/zone_grouping.R tests/testthat/test-zone-grouping.R
git commit -m "feat: add surveyzones_group_zones with starting_points anchoring"
```

---

### Task 9: End-to-end integration test, docs, and full suite

**Files:**
- Modify: `tests/testthat/test-zone-grouping.R` (append)
- Regenerate: `NAMESPACE`, `man/surveyzones_group_zones.Rd`, `man/surveyzones_zone_groups.Rd` (or similar, via roxygen2 — exact filenames depend on roxygen's title-derived naming)

**Interfaces:** none new — this task only adds test coverage and regenerates package metadata.

- [ ] **Step 1: Append the end-to-end integration test**

```r
test_that("surveyzones_group_zones end-to-end: shapes, D_max bound, labels, tract coverage", {
  skip_if_not_installed("ROI.plugin.glpk")

  # 7 zones across 2 partitions, positions on a line within each partition.
  # Partition A: zones at 0,1,2,10,11 (5 zones); D_max=5 groups {0,1,2} and {10,11}.
  # Partition B: zones at 0,50 (2 zones); D_max=5 -> 2 singleton groups.
  zone_positions <- c(
    A_Z1 = 0, A_Z2 = 1, A_Z3 = 2, A_Z4 = 10, A_Z5 = 11,
    B_Z1 = 0, B_Z2 = 50
  )
  zone_partition <- c(rep("A", 5), rep("B", 2))
  names(zone_partition) <- names(zone_positions)

  # 2 tracts per zone, both located at the zone's position
  za <- tibble::tibble(
    zone_id = rep(names(zone_positions), each = 2),
    tract_id = paste0(rep(names(zone_positions), each = 2), "_t", 1:2),
    partition_id = rep(zone_partition, each = 2)
  )
  tract_pos <- stats::setNames(rep(zone_positions, each = 2), za$tract_id)

  ids <- za$tract_id
  dists <- tidyr::expand_grid(origin_id = ids, destination_id = ids) |>
    dplyr::filter(origin_id != destination_id) |>
    dplyr::mutate(distance = abs(tract_pos[origin_id] - tract_pos[destination_id]))

  result <- surveyzones_group_zones(
    zone_assignments = za,
    sparse_distances = dists,
    D_max = 5,
    solver = "glpk",
    use_cache = FALSE,
    complete_distances = FALSE
  )

  expect_s3_class(result, "surveyzones_zone_groups")
  expect_true(all(c("zone_id", "partition_id", "center_tract_id", "n_tracts",
                     "group_id", "ggg", "zzz", "label") %in% names(result$zones)))
  expect_true(all(c("group_id", "partition_id", "center_zone_id", "n_zones",
                     "diameter", "ggg") %in% names(result$groups)))
  expect_true(all(c("zone_id", "tract_id", "visit_order") %in% names(result$sequence)))

  expect_equal(nrow(result$zones), 7L)
  expect_equal(sort(result$sequence$tract_id), sort(za$tract_id))

  # Label format
  expect_true(all(grepl("^\\d{3}\\.\\d{3}$", result$zones$label)))

  # ggg and zzz restart at 001 in each partition
  for (pid in c("A", "B")) {
    part_zones <- result$zones[result$zones$partition_id == pid, ]
    expect_true("001" %in% part_zones$zzz)
    expect_true("001" %in% part_zones$ggg)
  }

  # Partition A groups into 2 clusters (close triple + close pair)
  a_groups <- unique(result$zones$group_id[result$zones$partition_id == "A"])
  expect_equal(length(a_groups), 2L)

  # Partition B: two zones far apart -> 2 singleton groups
  b_groups <- unique(result$zones$group_id[result$zones$partition_id == "B"])
  expect_equal(length(b_groups), 2L)

  # No group spans two partitions: every zone in a given group_id has the
  # same partition_id
  partitions_per_group <- result$zones |>
    dplyr::summarise(n_partitions = dplyr::n_distinct(.data$partition_id), .by = "group_id")
  expect_true(all(partitions_per_group$n_partitions == 1L))
})

test_that("surveyzones_group_zones is deterministic across repeat calls", {
  skip_if_not_installed("ROI.plugin.glpk")
  za <- tibble::tibble(
    zone_id = paste0("Z", 1:6),
    tract_id = paste0("Z", 1:6),
    partition_id = c(rep("A", 4), rep("B", 2))
  )
  positions <- c(Z1 = 0, Z2 = 1, Z3 = 2, Z4 = 3, Z5 = 0, Z6 = 20)
  ids <- names(positions)
  dists <- tidyr::expand_grid(origin_id = ids, destination_id = ids) |>
    dplyr::filter(origin_id != destination_id) |>
    dplyr::mutate(distance = abs(positions[origin_id] - positions[destination_id]))

  args <- list(
    zone_assignments = za, sparse_distances = dists, D_max = 5,
    solver = "glpk", use_cache = FALSE, complete_distances = FALSE
  )
  first <- do.call(surveyzones_group_zones, args)
  second <- do.call(surveyzones_group_zones, args)

  expect_equal(first$zones, second$zones)
  expect_equal(first$groups, second$groups)
  expect_equal(first$sequence, second$sequence)
})
```

- [ ] **Step 2: Run the new tests to confirm they pass**

```bash
Rscript -e "devtools::load_all(quiet = TRUE); testthat::test_file('tests/testthat/test-zone-grouping.R')"
```

Expected: `[ FAIL 0 | ... ]` — 0 failures, across every test in the file.

- [ ] **Step 3: Regenerate documentation and NAMESPACE**

```bash
Rscript -e "devtools::document()"
```

Expected: creates `man/surveyzones_group_zones.Rd` and adds `export(surveyzones_group_zones)` plus `S3method(print,surveyzones_zone_groups)` to `NAMESPACE`. Review the diff to confirm no unrelated entries changed.

- [ ] **Step 4: Run the full package test suite**

```bash
Rscript -e "devtools::test()" 2>&1 | tail -30
```

Expected: `0 failures` across the whole suite (all pre-existing tests plus the new `test-zone-grouping.R` file).

- [ ] **Step 5: Commit**

```bash
git add tests/testthat/test-zone-grouping.R NAMESPACE man/
git commit -m "test: add end-to-end and determinism tests for surveyzones_group_zones; regenerate docs"
```

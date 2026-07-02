# Zone Grouping and Labeling Design

**Date:** 2026-07-02
**Status:** Approved

## Problem

`surveyzones_build_zones()` produces zones by solving a p-median problem over
tracts. But zones sometimes come from elsewhere — an external tract-to-zone
assignment, not built by this package's solver. Those zones have no computed
center, no travel-distance grouping, and no stable label.

We need a new entry point that, given an already-defined (possibly
partitioned) zone assignment:

1. Finds a center for each zone.
2. Groups zones into travel groups where every zone is within `D_max` of its
   group's center.
3. Assigns each zone a stable `ggg.zzz` label.
4. Produces a tract-level visit sequence, oriented across zone boundaries.

## New function: `surveyzones_group_zones()`

### Inputs

- `zone_assignments` — tibble/data.frame with `tract_id`, `zone_id`, optional
  `partition_id`. Every tract in a given `zone_id` must share the same
  `partition_id`; a zone spanning partitions is an error.
- `sparse_distances` — tract-level sparse distance tibble (`origin_id`,
  `destination_id`, `distance`), same format used throughout the package.
- `D_max` — numeric scalar. Distance bound for grouping: every zone in a
  group must be within `D_max` of that group's center zone.
- `solver` — character scalar, one of the existing supported solvers
  (`"glpk"` default, `"highs"`, `"cbc"`, `"symphony"`, `"spopt"`). Reused for
  the zone-grouping p-median.
- `max_time`, `rel_tol`, `verbose`, `use_cache` — pass-through to the reused
  p-median machinery, same defaults as `surveyzones_build_zones()`.
- `access_points`, `speed_kmh`, `complete_distances` — same distance
  completion pattern as `surveyzones_sequence()`: when
  `complete_distances = TRUE` (default), missing pairs are filled via
  haversine using `surveyzones_complete_distances()`. Required for both
  medoid-finding (needs complete within-zone distances) and tract
  orientation.
- `method`, `control` — seriation method/control passed to
  `surveyzones_sequence_tracts()` for tract-level ordering. Zone- and
  group-level TSPs always use `"TSP"` (matching `.sequence_zones()`'s
  existing behavior — not exposed as a parameter).
- `starting_points` — optional sf POINT object with a `partition_id`
  column, one row per partition. When a partition has a matching row, its
  `zzz` and `ggg` orders are anchored so the first zone and first group are
  the ones closest to that point (see Step 5.5 below). Partitions absent
  from `starting_points` keep the unanchored TSP orientation. Requires
  `access_points` regardless of `complete_distances`, since anchoring needs
  real coordinates — `sparse_distances` alone can't answer "distance from
  an arbitrary point to a tract."

No `max_workload_per_zone`-style cap is exposed: grouping is bounded by
`D_max` only, mirroring the uncapacitated default in `surveyzones_build_zones()`.

### Algorithm

**Step 1 — Partition.** Zones are grouped by `partition_id` (defaulting all
zones to a single `"all"` partition when the column is absent, consistent
with `surveyzones_partition()`). Everything below runs independently per
partition; groups and the two TSP passes never cross partition boundaries.

**Step 2 — Find zone centers (medoid).** For each zone, pick the tract that
minimizes the sum of distances to every other tract in the same zone
(unweighted — no workload column is required or used). Distances are
symmetrized by averaging `d(a,b)` and `d(b,a)`, matching `.seriate_1d()`.
Zones with a single tract trivially center on themselves. If any zone's
within-zone tract-pair distances are incomplete after optional haversine
completion, abort and list the affected zones — an incomplete medoid
candidate set can silently pick the wrong center, so this fails loud rather
than guessing.

**Step 3 — Zone-to-zone distance matrix.** Reuse the existing
`.build_zone_dist_matrix()` unchanged: `d(zone_A, zone_B) = min` pairwise
tract distance between members of A and B. This works directly against
`zone_assignments` (`tract_id`, `zone_id`) with no modification.

**Step 4 — Group zones within `D_max` (p-median, reused).** Treat zones as
the new atomic units and reuse the *existing* tract-level p-median solver
one level up:

- Build a synthetic "tracts" table: one row per zone in the partition,
  `tract_id = zone_id`, `expected_service_time = 1` (dummy weight, since
  grouping is uncapacitated).
- Melt the Step 3 distance matrix into the sparse long format
  (`origin_id`, `destination_id`, `distance`), dropping the diagonal and
  `Inf` entries.
- Call `surveyzones_build_zones_single()` with `D_max` = the grouping bound,
  `max_workload_per_zone = Inf`, and the requested `solver`.

Because pairs beyond `D_max` are excluded from the candidate set before the
MILP is built (the same mechanism that already bounds tract-to-center
distance in `surveyzones_build_zones()`), every zone lands within `D_max` of
its group's center zone by construction — no extra constraint needed.
Self-pairs guarantee feasibility even for a zone that is isolated beyond
`D_max` from everything else (it becomes a singleton group).

The solver's output is relabeled for this context: its "zone_id" (the
group's representative) becomes our `group_id`; its `center_tract_id`
becomes `center_zone_id`; its `n_tracts` becomes `n_zones`; its `diameter`
(now a zone-to-zone quantity) is kept as the group's diameter.

**Step 5 — Two independent partition-wide TSPs.** These do *not* need to be
consistent with each other — a group's members may land anywhere,
non-contiguously, in the zone TSP order.

- **Group TSP**: seriate the group centers (`center_zone_id` values) within
  the partition, using the Step 3 distance matrix restricted to those
  zones → group visit order → `ggg = sprintf("%03d", position)`.
- **Zone TSP**: seriate *all* zones in the partition (independent of group
  membership), using the full Step 3 distance matrix → `zzz =
  sprintf("%03d", position)`.

**Step 5.5 — Anchor to the partition's starting point (optional).** A TSP
path has exactly two valid global orientations (forward and reversed). When
`starting_points` supplies a row for this partition, both the zone order and
the group order are independently reoriented so each starts at whichever end
is closer to the starting point: compare the haversine distance from the
starting point to the *first* element's representative tract versus the
*last* element's, and reverse the whole sequence if the last is closer. This
picks the better of the two path endpoints without re-solving the TSP with a
fixed start node — it is not guaranteed to be the single globally-closest
zone if that zone happens to fall in the interior of the path, but for a
starting point outside/adjacent to the zone layout (the common case — a
field office near but not inside the surveyed area), the nearest zone is
typically a path endpoint anyway. Applied identically and independently to
both the zone order (zzz) and the group order (ggg).

**Step 6 — Labels.** For each zone: `label = paste0(ggg, ".", zzz)`, where
`ggg` comes from the zone's group (Step 5's group TSP) and `zzz` from the
zone itself (Step 5's zone TSP). Labels are unique within a partition, not
globally (no partition prefix, per requirements) — callers needing global
uniqueness keep `partition_id` as a separate column.

**Step 7 — Tract-level sequencing and cross-zone orientation.** The zone TSP
order from Step 5 doubles as the walk path for tract orientation: process
zones in `zzz` order (within each partition), and for each zone:

1. Order its tracts with `surveyzones_sequence_tracts()` (reused unchanged —
   it already works on a bare tract-ID vector plus `sparse_distances`).
2. Orient the resulting list (reverse if needed) so the zone's entry tract
   is closest to the previous zone's exit tract, walking in `zzz` order.
   This is the same mechanic as `.sequence_tracts_with_orientation()`,
   adapted to walk the zone-TSP order instead of a group+within-group order.

### Output: new S3 class `surveyzones_zone_groups`

- `zones` — `zone_id` (original input ID), `partition_id`, `center_tract_id`,
  `n_tracts`, `group_id`, `ggg`, `zzz`, `label`
- `groups` — `group_id`, `partition_id`, `center_zone_id`, `n_zones`,
  `diameter`, `ggg`
- `sequence` — `zone_id`, `tract_id`, `visit_order` (tract-level, oriented,
  following the `zzz` walk order)
- `parameters` — `D_max`, `solver`, `timestamp`
- `diagnostics` — solver status/objective/solve_time per partition, from
  Step 4's reused solve

### Error handling

- `zone_assignments` missing required columns, or a `tract_id` appearing in
  more than one `zone_id` → abort (mirrors `validate_tracts()`-style checks
  elsewhere).
- A zone whose tracts span more than one `partition_id` → abort, listing the
  offending zones.
- Incomplete within-zone distances for medoid-finding (Step 2), after
  optional haversine completion → abort, listing affected zones.
- Unsupported `solver` → reuse `validate_solver()` unchanged.
- `starting_points` supplied without `access_points` → abort (anchoring
  needs real coordinates).
- `starting_points` supplied without a `partition_id` column → abort.

### Testing

Organized by the step each test targets, so the implementation plan can
drive them test-first per step.

**Center-finding (Step 2 — medoid)**
- 3-tract zone with hand-computed distances: assert the medoid is the tract
  minimizing total distance to the other two.
- Single-tract zone: center is that tract, no distance lookup needed.
- Asymmetric input (`d(a,b) != d(b,a)`): assert the averaged-symmetrized
  value is what's used, via a case where the naive (non-averaged) choice
  would pick a different tract.
- Exact tie in total distance between two candidate tracts: assert
  deterministic selection (documented tie-break, e.g. lowest `tract_id`),
  not solver/platform-dependent.
- Incomplete within-zone distances (a missing pair) with
  `complete_distances = FALSE`: assert abort naming the specific zone.
- Two different zones both missing pairs: assert the abort message lists
  *both* zone IDs, not just the first encountered.
- `complete_distances = TRUE` with `access_points` supplied: a zone missing
  a pair in `sparse_distances` still resolves a medoid, using the
  haversine-filled value.

**Starting-point anchoring (Step 5.5)**
- Given a starting point clearly closer to one end of a known line layout:
  assert the zone whose representative tract is nearest the starting point
  ends up `zzz == "001"`, and its group ends up `ggg == "001"`.
- Given no `starting_points` row for a partition: assert that partition's
  order is unaffected (same as the no-anchoring case).
- Given a single zone or single group: anchoring is a no-op (nothing to
  reverse).
- `starting_points` provided but missing `access_points` → abort.

**Zone grouping (Step 4 — reused p-median)**
- Property test: for every zone in the output, `distance(zone, group's
  center_zone_id) <= D_max` — generate several random small synthetic
  layouts and assert this holds in all of them, not just one hand-built
  case.
- Multi-partition input: assert no group contains zones from two different
  `partition_id` values.
- Isolated zone (all pairwise distances to other zones exceed `D_max`):
  assert it becomes a singleton group (`n_zones == 1`), and the run
  succeeds rather than reporting infeasible.
- Two well-separated clusters of zones (small intra-cluster distances,
  large inter-cluster distances, `D_max` set between the two scales):
  assert exactly 2 groups form, one per cluster.
- Single-zone partition: assert it short-circuits to one trivial group
  without error.
- `solver = "spopt"` produces the same grouping structure as `"glpk"` on a
  small case (skip via `testthat::skip_if_not_installed("spopt")`).
- `use_cache = TRUE`: second call with identical arguments returns an
  identical result (mirrors the existing cache-hit tests for
  `surveyzones_build_zones()`).

**Labeling (Step 5–6 — two independent TSPs)**
- Every `label` matches `^\d{3}\.\d{3}$`.
- `zzz` values within a partition form a permutation of `1..n_zones` (no
  gaps, no duplicates).
- `ggg` values within a partition form a permutation of `1..n_groups`.
- Explicit non-contiguity case: construct a layout where a group's members
  are *not* adjacent in the zone TSP order, and assert the label output
  reflects that (e.g. group `002` contains zones with `zzz` values `003`,
  `007`, `014`) — this is a regression guard against accidentally coupling
  the two TSPs.
- Multi-partition input: assert both `ggg` and `zzz` restart at `001` in
  each partition independently.
- Labels may repeat across partitions (no global uniqueness claim) but must
  be unique within a partition.

**Tract sequencing and orientation (Step 7)**
- `sequence` output contains exactly the same set of `tract_id` values as
  the input `zone_assignments`, no more, no fewer.
- `visit_order` is a `1..n` permutation within each zone.
- Orientation regression case: two adjacent zones (in `zzz` order) where
  reversing the second zone's tract order strictly reduces the entry
  distance from the first zone's exit tract — assert the reversal happens
  (mirrors the intent of `.sequence_tracts_with_orientation()`'s existing
  coverage).
- Orientation does not carry across partition boundaries: the first zone of
  partition 2 is not oriented relative to the last zone of partition 1.
- Zones with 1–2 tracts (where `surveyzones_sequence_tracts()` short-circuits
  and returns tracts unchanged) still participate correctly in orientation.

**Validation**
- Missing `tract_id` or `zone_id` column: informative abort naming the
  missing column.
- A `tract_id` assigned to two different `zone_id` values: abort.
- A zone whose tracts span two different `partition_id` values: abort,
  naming the zone.
- Unsupported `solver` value: abort via `validate_solver()`, message
  matches existing behavior.
- `complete_distances = TRUE` with no `access_points` available (neither
  passed nor inferable): abort with the same guidance message pattern as
  `surveyzones_sequence()`.

**End-to-end integration**
- One synthetic dataset (~10–15 tracts, 4–5 zones, 2 partitions) exercised
  through the full function: assert the shapes and column sets of all five
  output components (`zones`, `groups`, `sequence`, `parameters`,
  `diagnostics`) together, plus the cross-cutting invariants above (D_max
  bound, label format, tract coverage) all hold simultaneously — catches
  integration bugs that per-step unit tests could miss (e.g. a rename
  applied to the wrong table).
- Determinism: running the same inputs twice (same solver, same seed where
  relevant) yields identical `zones`/`groups`/`sequence` output — not just
  structurally-valid output that happens to differ run to run.

### What stays the same (reused, unmodified)

- `surveyzones_partition()`, `surveyzones_build_zones_single()`,
  `validate_solver()`, `.build_zone_dist_matrix()`, `.seriate_1d()`,
  `surveyzones_sequence_tracts()`, `surveyzones_complete_distances()`.
- `sf::st_distance()` — used directly for starting-point anchoring (Step
  5.5), same pattern `surveyzones_complete_distances()` already uses for
  haversine fill-in. No new package dependency (`sf` is already imported).

### Out of scope

- No workload/size cap on groups (distance-only, per requirements).
- No changes to `surveyzones_build_zones()`, `surveyzones_sequence()`, or the
  existing `{partition}_{group}.{pos:03d}` zone ID renaming — this is a new,
  separate entry point for externally-defined zones, not a replacement.
- No reconciliation between the group TSP and zone TSP orders — they are
  independent by design (starting-point anchoring is applied to each
  independently, and does not make them agree with each other).
- No exact "true nearest zone is first" guarantee — anchoring picks the
  better of the TSP path's two endpoints (see Step 5.5), not a full
  re-optimization with a fixed start node.

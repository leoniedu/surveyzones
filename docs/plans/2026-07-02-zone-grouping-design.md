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

### What stays the same (reused, unmodified)

- `surveyzones_partition()`, `surveyzones_build_zones_single()`,
  `validate_solver()`, `.build_zone_dist_matrix()`, `.seriate_1d()`,
  `surveyzones_sequence_tracts()`, `surveyzones_complete_distances()`.

### Out of scope

- No workload/size cap on groups (distance-only, per requirements).
- No changes to `surveyzones_build_zones()`, `surveyzones_sequence()`, or the
  existing `{partition}_{group}.{pos:03d}` zone ID renaming — this is a new,
  separate entry point for externally-defined zones, not a replacement.
- No reconciliation between the group TSP and zone TSP orders — they are
  independent by design.

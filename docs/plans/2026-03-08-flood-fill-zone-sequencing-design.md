# Flood-Fill Zone Sequencing Design

**Date**: 2026-03-08
**Status**: Approved

## Problem

Seriation-based zone sequencing (SPIN_NH, TSP, etc.) produces a global 1D
ordering that can interleave distant zones to optimize the overall path.
Fieldworkers need spatially coherent groups of consecutive zones — everything
close should be visited before jumping far.

## Solution

Replace seriation with a two-level **flood-fill + TSP** approach:

1. **Build zone distance matrix** — reuse `.build_zone_dist_matrix()` (min
   pairwise tract distances) on the completed distance table.
2. **Flood-fill to find groups** — from an unvisited zone, collect all
   unvisited zones within `threshold` (transitive closure). Produces groups
   of threshold-connected zones.
3. **Seriate within each group** — TSP seriation on zones within each group
   (via `.seriate_1d()` on the zone distance matrix subset). Groups of 1–2
   zones skip seriation.
4. **Seriate between groups** — build group-to-group distance matrix:
   `d(group_A, group_B) = min(d(zone_i, zone_j))` for all `i in A, j in B`.
   TSP-seriate the groups.
5. **Concatenate** — group order × within-group order → final `zone_order`.

## Interface

```r
surveyzones_sequence_zones(
  plan,
  sparse_distances,
  access_points = NULL,
  speed_kmh = 0.1,
  complete_distances = TRUE,
  threshold = 10,
  by_partition = TRUE
)
```

- **Removed**: `method`, `control` — TSP is used internally, not exposed.
- **Added**: `threshold` — numeric scalar in distance units (minutes).
  Default `10`.
- **Unchanged**: `access_points`, `speed_kmh`, `complete_distances`,
  `by_partition`.
- **Output**: same `plan$zone_sequence` tibble (`partition_id`, `zone_id`,
  `zone_order`).

## New internals

- `.flood_fill_groups(mat, threshold)` — takes a zone distance matrix and
  threshold, returns a list of character vectors (zone IDs per group).

## What stays the same

- `.build_zone_dist_matrix()`, `.seriate_1d()`, `surveyzones_sequence()`,
  `surveyzones_sequence_tracts()`, distance completion flow — all unchanged.
- `seriation` package still used internally (`.seriate_1d()` with `"TSP"`
  for within-group and between-group ordering).

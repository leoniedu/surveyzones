# Gap-Based Zone Grouping Design

**Date**: 2026-03-08
**Status**: Approved

## Problem

Flood-fill zone grouping uses transitive connectivity: if A is within
threshold of B, and B within threshold of C, then A-B-C form one group.
In dense urban areas, even with a small threshold (5 min), this creates
one giant group containing all zones — the `group_id` becomes meaningless.

## Solution

Replace flood-fill with **TSP + gap splitting**:

1. **Build zone distance matrix** — same `.build_zone_dist_matrix()` (min
   pairwise tract distances), unchanged.
2. **TSP-seriate all zones** — single `.seriate_1d()` call on all zones in
   the partition, producing one optimal spatial path.
3. **Scan consecutive distances** — walk the TSP path, look up
   `d(zone[i], zone[i+1])` from the zone distance matrix.
4. **Split at gaps** — whenever consecutive distance > `threshold`,
   increment `group_id`.

## Why this is better

- **No transitive blowup** — groups are defined by actual gaps in the TSP
  path, not transitive connectivity.
- **Simpler** — one TSP pass + linear scan, instead of flood-fill +
  within-group TSP + between-group TSP + concatenation.
- **Same zone_order** — the TSP path determines visit order regardless of
  grouping. `group_id` is a label on top.

## Interface

Same public API, same output schema:

```r
surveyzones_sequence_zones(
  plan, sparse_distances,
  access_points = NULL, speed_kmh = 0.1,
  complete_distances = TRUE, threshold = 10,
  by_partition = TRUE
)
```

- `threshold` changes meaning: "max consecutive distance before group
  break" (was: flood-fill connectivity radius).
- Output: `partition_id`, `zone_id`, `zone_order`, `group_id`.

## What changes

- `.flood_fill_groups()` — **removed**.
- Inner logic of `surveyzones_sequence_zones()` — simplified to single TSP
  + gap scan.

## What stays the same

- `.build_zone_dist_matrix()`, `.seriate_1d()`, `.matrix_to_sparse()` —
  unchanged.
- `surveyzones_sequence()`, `surveyzones_sequence_tracts()` — unchanged.
- Distance completion flow — unchanged.
- pns.zonas integration — unchanged (same API).

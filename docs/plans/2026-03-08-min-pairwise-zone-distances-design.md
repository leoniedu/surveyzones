# Min-Pairwise Zone Distance for Zone Sequencing

**Date**: 2026-03-08
**Status**: Approved

## Problem

`surveyzones_sequence_zones()` sequences zones using center-to-center distances.
Two zones whose centers are far apart may have border tracts that are practically
touching. This causes visually adjacent zones to receive distant sequence numbers.

## Solution

Replace center-to-center distances with **minimum pairwise tract distances** when
building the zone-to-zone distance matrix for seriation.

`d(zone_A, zone_B) = min(d(tract_i, tract_j))` for all `i in A, j in B`.

## Design

### New helper: `.build_zone_dist_matrix()`

Builds a symmetric zone-to-zone distance matrix from the sparse distance table:

1. Join `sparse_distances` with `plan$assignments` to tag each origin/destination
   with its `zone_id`.
2. Filter to cross-zone pairs (`zone_origin != zone_dest`).
3. Group by `(zone_origin, zone_dest)`, take `min(distance)`.
4. Symmetrise: `d(A,B) = min(d(A,B), d(B,A))`.
5. Return a named numeric matrix (zone_ids as row/col names).

### Changes to `surveyzones_sequence_zones()`

- Call `.build_zone_dist_matrix()` instead of passing center IDs to `.seriate_1d()`.
- Build the zone distance matrix directly, then call `seriation::seriate()` on it.
- For missing zone pairs (no cross-zone tract pair in sparse table), fall back to
  haversine between access points: compute `min(haversine(pt_i, pt_j))` for all
  `pt_i in zone_A, pt_j in zone_B` using the `access_points` sf object, penalised
  by `speed_kmh`.
- The `complete_distances`, `access_points`, and `speed_kmh` parameters keep the
  same interface but now fill missing zone-pair distances rather than center-pair
  distances.

### What stays the same

- `surveyzones_sequence()`, `surveyzones_sequence_tracts()`, `.seriate_1d()` unchanged.
- `method`, `control`, `by_partition` parameters work identically.
- Output schema (`plan$zone_sequence`) identical.

## Trade-offs

- Requires cross-zone tract pairs in the sparse distance table. These exist for
  nearby zones (within D_max). Distant zone pairs fall back to haversine.
- Slightly more computation than center-to-center (aggregate over all cross-zone
  pairs), but zone counts are small (tens, not thousands).

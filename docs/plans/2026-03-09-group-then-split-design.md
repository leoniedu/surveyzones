# Two-Phase Group-Then-Split Zone Building

*Design date: 2026-03-09*

## Problem

The current `uncap_then_split` strategy initialises the uncapacitated phase with
`K₀ = max(K_from_size, K_from_workload)`. This means workload considerations inflate
the number of geographic groups produced in phase 1, creating finer clusters than the
geography alone requires. The new algorithm lets geography speak first.

## Algorithm

### Phase 1 — geographic grouping (uncapacitated)

Run the existing K-search loop with `max_workload_per_zone = Inf`, starting at **K=1**.
Increment K until all tracts are assignable within D_max. The resulting zone assignments
become the `group_id`. This is a purely geographic operation — workload plays no role.

### Phase 2 — workload splitting (capacitated, per group)

For each phase-1 group:

- If `total_workload ≤ max_workload_per_zone`: keep as a single zone, assign `group_id`.
- If `total_workload > max_workload_per_zone`: run a fresh K-search on that group's tracts
  starting at **K=1**, with both D_max and workload constraints active. All resulting
  sub-zones inherit the parent `group_id`.

Phase 2 sub-problems start at K=1 (not `ceiling(wl / max_workload)` as in the old code)
and search upward, following the same pattern as phase 1.

## API Changes

### `target_zone_size` removed

The parameter was used to inflate K₀ and is no longer needed. Removed from:
- `surveyzones_build_zones()`
- `surveyzones_build_zones_single()`
- Internal recursive calls

### `strategy` values

| Value | Behaviour |
|-------|-----------|
| `"auto"` + capacitated | Two-phase group-then-split (new default) |
| `"auto"` + uncapacitated | Single-phase K-search from K=1 |
| `"direct"` | Single MILP, K₀ = ceiling(total_workload / max_workload_per_zone) |

`"uncap_then_split"` is removed as a valid strategy value.

## Data Model

### New `group_id` column

Added to `plan$assignments` and `plan$zones` (character):

- Identifies which phase-1 geographic group a tract/zone belongs to.
- For unsplit groups: `group_id = zone_id`.
- For split groups: all child zones share the parent `group_id`.
- Uncapacitated case: `group_id = zone_id` (each zone is its own group).

### Diagnostics

Add `n_groups` (integer) to `plan$diagnostics` — number of phase-1 groups.

## Error Handling

- Phase 1 infeasible → fall back to `"direct"` with a warning (same as today).
- Phase 2 sub-problem infeasible → keep oversized zone, emit warning (same as today).

## Downstream Impact

- `surveyzones_sequence()`: the existing `group_id` column in `plan$zone_sequence`
  (used for gap-based grouping) is separate from the new zone-level `group_id`.
  No sequencing changes required.
- `pns.zonas`: callers using `target_zone_size` will need to remove that argument.
  The `strategy` argument can be removed from call sites since `"auto"` now does the
  right thing.

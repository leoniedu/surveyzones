# Sequential Zone ID Design

**Date:** 2026-03-09
**Status:** Approved

## Problem

Current `zone_id` values are derived from phase-1 center tract IDs (e.g.
`280670105000109.001`). These are:

- **Not anonymous** — embed internal tract identifiers
- **Not sequential** — bear no relation to visit order
- **Conflicting** — `assignments$group_id` (character, phase-1 center) collides
  with `zone_sequence$group_id` (integer, travel group) in downstream joins

## Design

### ID Format

```
{partition_id}_{group_id}.{pos:03d}
```

| Component    | Description                                                      | Example |
|--------------|------------------------------------------------------------------|---------|
| `partition_id` | Existing single-column partition value (e.g. UF code)          | `28`    |
| `group_id`   | Integer travel group within that partition (1, 2, 3…)           | `2`     |
| `pos`        | 1-indexed position within `(partition_id, group_id)`, 3-digit   | `003`   |

Examples: `28_1.001`, `28_2.003`, `11_3.007`

IDs are **globally unique** (partition prefix disambiguates across partitions).
Within a partition they sort to reflect travel sequence order.

### Implementation Location

Applied at the end of `surveyzones_sequence()`, after both `zone_sequence` and
`sequence` are computed — the only point where all required information
(partition, travel group, zone order) is available together.

### Rename Step

A new internal helper `.rename_zones(plan)`:

1. Build lookup from `zone_sequence`:
   - `pos = row_number()` within each `(partition_id, group_id)`, ordered by
     `zone_order`
   - `new_zone_id = paste0(partition_id, "_", group_id, ".", sprintf("%03d", pos))`
2. Apply old → new mapping to:
   - `plan$assignments$zone_id`
   - `plan$zones$zone_id`
   - `plan$zone_sequence$zone_id`
   - `plan$sequence$zone_id`
3. Drop `plan$assignments$group_id` (character, phase-1 center — redundant with
   `plan$zones$center_tract_id` which is preserved unchanged)

### Downstream fix: `pns.zonas/R/map.R`

The `.travel_group` workaround introduced to avoid the `group_id` naming
conflict is removed. With `assignments$group_id` dropped, the join is clean and
`color_by(.data$group_id)` can be used directly.

## What Is Preserved

- `plan$zones$center_tract_id` — original center tract, unchanged
- `plan$zones$partition_id`, `plan$assignments$partition_id` — unchanged
- All other columns — unchanged

## Out of Scope

- Renaming zones when `surveyzones_sequence()` is not called (plan stays with
  old IDs; caller must sequence before using final IDs)
- Supporting alternative ID formats (can be added later as a parameter)

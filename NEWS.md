# surveyzones 0.0.0.9003

* `surveyzones_assign_teams(use_cache = TRUE)` memoises the whole solve on disk
  (separate `assign_teams` cache dir, `verbose` not in the key). On a hit,
  `diagnostics$solve_time` is that of the original solve.
  `surveyzones_clear_cache()` clears it and now keeps the cache directories
  (removing them broke later memoised calls in the same session).
  After upgrading surveyzones, run `surveyzones_clear_cache()`: cache keys do
  not include helper code.

# surveyzones 0.0.0.9002

* `surveyzones_assign_teams(k = "auto")` picks `k` from the minimum number of
  roteiros (`k_min = ceiling(min_roteiros / n_teams)`, `k = k_min + 2`);
  `parameters` records `k`, `k_auto` and `k_min`.
  Non-whole `k` is now rejected.
* Solver status handling recognises glpk status symbols (`GLP_OPT`,
  `GLP_NOFEAS`, `GLP_INFEAS`, ...), fixing optimal glpk solves being reported
  as merely feasible.

# surveyzones 0.0.0.9001

* New `surveyzones_roteiro_candidates()` builds candidate walkable roteiros
  (sets of tracts reachable on foot from a centre household) for team
  assignment.
* New `surveyzones_assign_teams()` assigns every tract to one team with a
  sparse mixed-integer model that chooses roteiros, loose tracts and team
  origins under load-balance, roteiro-count and coverage rules. Reports an
  elastic diagnosis and the minimum number of roteiros when infeasible.

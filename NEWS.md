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

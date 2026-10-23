# surveyzones 0.0.0.9001

* New `surveyzones_roteiro_candidates()` builds candidate walkable roteiros
  (sets of tracts reachable on foot from a centre household) for team
  assignment.
* New `surveyzones_assign_teams()` assigns every tract to one team with a
  sparse mixed-integer model that chooses roteiros, loose tracts and team
  origins under load-balance, roteiro-count and coverage rules. Reports an
  elastic diagnosis and the minimum number of roteiros when infeasible.

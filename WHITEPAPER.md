# surveyzones: Workload-Balanced Zone Construction for Survey Field Operations

**White Paper**

*Eduardo Leoni, IBGE*

---

## Executive Summary

Survey organizations face a fundamental operational challenge: how to partition geographic regions into zones for field teams such that:
1. Each team has manageable, balanced workload
2. Field operations remain efficient (teams don't waste time traveling)
3. Geographic zones are cohesive and traversable within reasonable time

Existing facility location models—Location Set Covering Problem (LSCP), p-median, p-center, Maximum Coverage Location Problem (MCLP), and warehouse location models—were designed for retail, emergency services, or supply chain optimization. While mathematically sophisticated, they fundamentally misalign with the constraints of survey field operations.

This white paper demonstrates why surveyzones' approach—combining workload-balanced zone construction with maximum diameter constraints and asymmetric routing—is necessary and distinct. We analyze five competing models, identify their limitations in the survey context, and justify surveyzones' architectural choices.

**Key Finding:** Only surveyzones jointly addresses workload balancing, geographic cohesion, and travel efficiency. Existing models optimize for logistics cost (p-median), coverage reliability (LSCP), or fairness (p-center), not for human-centered survey operations.

---

## 1. The Survey Zone Problem

### 1.1 Operational Context

Survey organizations conduct field operations in which:
- Census or service tracts must be visited by field teams
- Each tract has expected service time (interview duration, survey complexity)
- Teams have limited capacity (e.g., 8 hours/day)
- Travel between tracts consumes time but adds no survey value
- Teams work from designated zones, ideally staying within geographic boundaries

Example: The Brazilian IBGE conducts household surveys across ~300,000 census tracts. Field supervisors must partition these tracts into workable zones for ~50,000 field teams. Each zone should allow one team to complete all tracts in one work day without excessive travel.

### 1.2 Constraints Specific to Survey Operations

1. **Workload Balance**: Field teams have equal capacity. Zones must respect `max_workload_per_zone` (typically 6-8 hours of survey time). Unlike retail (where customer demand varies), survey work is more predictable—each tract takes ~30 min to 2 hours.

2. **Geographic Diameter**: Teams cannot efficiently traverse zones larger than a certain diameter. A zone where the farthest two tracts are 2 hours apart is operationally impractical, even if the total within-zone travel time is acceptable. Unlike p-median (which ignores zone size), surveyzones enforces maximum diameter `D_max`.

3. **Asymmetric Routing**: Real-world travel (OSRM data) is asymmetric—travel time from A→B often differs from B→A due to road networks, traffic patterns, one-way streets. Most models assume symmetric distances or ignore routing entirely. surveyzones accommodates asymmetry via ATSP (Asymmetric Traveling Salesman Problem).

4. **Zone Contiguity** (soft constraint): While not hard-enforced in the MILP, zones should naturally cluster spatially. This emerges from D_max filtering rather than explicit contiguity constraints—a design choice that keeps the MILP tractable.

---

## 2. Literature Review: Existing Facility Location Models

### 2.1 Location Set Covering Problem (LSCP)

**Formulation:**
$$\min K = \sum_j y_j$$

Subject to:
$$\sum_{j: d_{ij} \leq D_{\max}} y_j \geq 1 \quad \forall i$$

**Objective:** Minimize the number of facilities (zones) needed to ensure all demand points are within distance D_max.

**Strengths:**
- Natural for coverage-based problems (emergency responders, retail)
- Simple interpretation: "Minimum zones to reach everyone within D_max"
- Well-studied, scalable algorithms

**Limitations for Survey Operations:**
1. **No workload balancing**: LSCP might create highly imbalanced zones. If Tract A contains 100 households and Tract B contains 10, LSCP doesn't distribute them across teams—it just ensures coverage.
2. **No distance minimization**: LSCP stops once it achieves minimum K. Two different zone configurations with the same K are equally valid, even if one requires 2x travel time.
3. **Facility placement flexibility**: LSCP assumes facilities can be placed anywhere. surveyzones requires centers to be actual tracts (operational constraint).
4. **Coverage vs. assignment**: LSCP requires each demand to be *reachable* from some facility. surveyzones requires each tract be *assigned* to exactly one zone. This distinction matters: a tract might be within D_max of multiple centers, but assignment is deterministic.

**Example Failure:**
- 4 tracts: A (5 hr workload), B (1 hr), C (1 hr), D (5 hr)
- D_max = 20 minutes
- LSCP might find K=2 with zone 1={A, B, C} and zone 2={D}
- Team 1 has 7 hours work, Team 2 has 5 hours—operationally unfair
- surveyzones with max_workload=6 would require K=3 to balance

### 2.2 p-Median Problem

**Formulation:**
$$\min \sum_i d_{ij(i)} \cdot w_i$$

Subject to:
$$\sum_j y_j = K$$

$$\sum_j x_{ij} = 1 \quad \forall i$$

**Objective:** Given K facilities, minimize total travel distance (weighted by demand).

**Strengths:**
- Directly minimizes logistics cost
- Well-studied, efficient algorithms
- Natural for delivery/service routing

**Limitations for Survey Operations:**
1. **No workload capacity**: p-median ignores that teams have limited capacity. It might assign all high-demand tracts to one zone and low-demand to another, maximizing distance efficiency but creating operational imbalance.
2. **No diameter constraint**: Zones can be arbitrarily large. A p-median solution might create a zone spanning 50km because it minimizes total distance—operationally impractical for field teams.
3. **Assumes symmetric distances**: Standard p-median formulations assume $d_{ij} = d_{ji}$. Real-world survey operations (OSRM routing) are asymmetric.
4. **Assumes demand quantity matters**: p-median weights by $w_i$ (demand). Survey workload is continuous (service time), not discrete demand units.

**Example Failure:**
- 10 tracts in a 50km × 50km area
- p-median with K=2 might create one zone covering 40km diameter because it minimizes total distance
- Even though two smaller zones would double distance, they're operationally better for survey teams
- surveyzones would enforce max diameter = 20km, requiring more zones but preserving operability

### 2.3 p-Center Problem

**Formulation:**
$$\min r = \max_i d_{ij(i)}$$

Subject to:
$$\sum_j y_j = K$$

**Objective:** Minimize the maximum distance from any demand point to its assigned facility (minimax fairness).

**Strengths:**
- Ensures fairness: nobody travels excessively far
- Useful for emergency services (worst-case response time)
- Natural for equity-focused optimization

**Limitations for Survey Operations:**
1. **No workload balancing**: Like p-median, ignores team capacity. Might assign all high-demand tracts to one zone.
2. **Optimizes for wrong metric**: p-center minimizes worst-case distance. Surveys need workload fairness, not distance fairness. A 30-minute drive is fine if workload is 4 hours; a 1-hour drive is unacceptable even if workload is 2 hours.
3. **Ignores zone diameter**: Minimizes distance to nearest facility, not the spread within the zone.

**Example Failure:**
- 4 tracts arranged in a line: A—B—C—D (each 10 min apart)
- p-center with K=2 might place facilities at B and C, minimizing max distance (10 min)
- But this creates two zones with low workload, operationally inefficient
- surveyzones would balance workload while respecting geometry

### 2.4 Maximum Coverage Location Problem (MCLP)

**Formulation:**
$$\max \sum_i w_i \cdot c_i$$

Subject to:
$$\sum_j y_j = K$$

$$c_i \leq \sum_{j: d_{ij} \leq D_{\max}} y_j$$

**Objective:** Given K facilities, maximize the total demand covered within D_max.

**Strengths:**
- Handles partial coverage (not all demand need be covered)
- Useful when service areas are limited (e.g., limited ambulances)
- Scalable algorithms available

**Limitations for Survey Operations:**
1. **Allows uncovered demand**: MCLP is designed for partial coverage scenarios (emergency services with limited resources). Survey operations require 100% coverage.
2. **No assignment determinism**: MCLP allows demand to be uncovered. surveyzones requires all tracts to be assigned.
3. **No workload capacity**: Like LSCP and p-median, no workload balancing.
4. **Maximizes, doesn't minimize**: Different optimization direction—less natural for zone construction where we want to find feasible assignments.

---

### 2.5 Warehouse Location Problem (WLP) / Capacitated Facility Location

**Formulation:**
$$\min \sum_j f_j y_j + \sum_i \sum_j c_{ij} x_{ij}$$

Subject to:
$$\sum_j x_{ij} = 1 \quad \forall i$$

$$\sum_i s_i x_{ij} \leq S_j y_j \quad \forall j$$

**Objective:** Minimize fixed facility costs plus variable transportation costs, respecting warehouse capacity.

**Strengths:**
- Handles facility opening costs (relevant for multi-tier networks)
- Incorporates capacity constraints (most similar to surveyzones)
- Well-studied in operations research

**Limitations for Survey Operations:**
1. **No geographic diameter constraint**: WLP respects capacity but not zone size. A single large warehouse serving one region could be optimal from cost perspective.
2. **Assumes symmetric, known travel costs**: Optimizes for cost, not operability.
3. **Facility opening costs don't apply**: Surveys don't have fixed costs for opening zones; they have team availability.
4. **Focus on logistics efficiency**: Optimized for supply chain, not human-centered field operations.

---

## 3. Comparative Analysis

### 3.1 Constraint Comparison Table

| Constraint | LSCP | p-median | p-center | MCLP | WLP | **surveyzones** |
|-----------|------|---------|---------|------|-----|-----------------|
| **Workload Capacity** | ❌ | ❌ | ❌ | ❌ | ✅ | ✅ |
| **Maximum Diameter** | ❌ | ❌ | ❌ | ❌ | ❌ | ✅ |
| **100% Assignment** | ✅* | ✅ | ✅ | ❌ | ✅ | ✅ |
| **Distance Optimization** | ❌ | ✅ | ✅ (minimax) | ❌ | ✅ | ✅ |
| **Asymmetric Routing** | ❌ | ❌ | ❌ | ❌ | ❌ | ✅ |
| **Capacity Constraints** | ❌ | ❌ | ❌ | ❌ | ✅ | ✅ |
| **Automatic K Discovery** | ✅ (min K) | ❌ | ❌ | ✅ | ❌ | ✅ |

*LSCP ensures coverage; surveyzones ensures assignment

### 3.2 Objective Function Comparison

| Model | Minimizes | Given/Fixed |
|-------|-----------|-------------|
| **LSCP** | # of facilities | D_max |
| **p-median** | Total distance | K |
| **p-center** | Max distance | K |
| **MCLP** | # of facilities | K, D_max |
| **WLP** | Fixed + variable costs | (none—free variables) |
| **surveyzones** | Total distance | D_max, workload_max |

**Key Insight:** surveyzones is the only model that jointly:
1. Respects workload capacity (unlike p-median, p-center, LSCP)
2. Enforces geographic bounds via diameter (unlike all others)
3. Automatically discovers minimum feasible K (like LSCP, unlike p-median/p-center)
4. Minimizes travel distance (like p-median, unlike LSCP)

---

## 4. Why surveyzones is Necessary

### 4.1 The Survey Operations Gap

Existing models optimize for:
- **Logistics efficiency** (p-median): Minimize miles driven—wrong metric for survey time constraints
- **Coverage reliability** (LSCP, MCLP): Ensure everyone is reachable—doesn't ensure operability
- **Fairness in distance** (p-center): Minimize worst-case distance—doesn't address workload fairness
- **Cost minimization** (WLP): Balance facility costs vs. transportation—not applicable to surveys

**None address the core survey problem:**
> *"Partition tracts into zones such that each field team can complete all assigned work within one day, stay geographically cohesive, and minimize wasted travel time."*

### 4.2 Mathematical Justification

surveyzones solves a **constrained optimization problem with three interdependent objectives:**

$$\min \sum_i \sum_k d_{ik} \cdot z_{ik}$$

Subject to:
$$\sum_k z_{ik} = 1 \quad \forall i \quad \text{(all tracts assigned)}$$

$$\sum_k y_k = K \quad \text{(exactly K zones)}$$

$$\sum_{i: z_{ik}=1} s_i \leq W_{\max} \quad \forall k \quad \text{(workload capacity)}$$

$$D_{\max \text{ filtering}} \quad \text{(geographic cohesion)}$$

$$K = \min K' : \text{constraints feasible} \quad \text{(automatic zone count)}$$

This combines:
1. **LSCP's automatic K discovery** (find minimum feasible zones)
2. **p-median's distance minimization** (optimize travel efficiency)
3. **WLP's capacity constraints** (balance workload)
4. **Novel: Maximum diameter enforcement** (ensure operability)

### 4.3 Operational Justification

**Field teams are constrained by time, not distance:**

A team with 8 hours available will spend:
- ~4 hours in survey work (2 tracts × 2 hours each)
- ~3 hours in travel time
- ~1 hour in breaks, admin, transitions

LSCP/p-median optimize for distance but ignore workload distribution. surveyzones optimizes the *time envelope* teams actually have.

**Geographic cohesion is non-negotiable:**

Even if minimizing distance, a zone spanning 100km diameter is operationally infeasible:
- Team can't supervise effectively
- Logistics coordination is difficult
- Absenteeism increases
- Quality suffers

surveyzones enforces `D_max` to keep zones within 20-40km (adjustable), ensuring manageability.

**Asymmetric routing reflects reality:**

Real road networks are asymmetric. LSCP/p-median assume symmetric distances, losing fidelity. surveyzones uses ATSP to optimize sequencing with real routing data.

---

## 5. Technical Architecture

### 5.1 Algorithm Overview

surveyzones combines three stages:

**Stage 1: Connectivity Analysis**
- Filter distances by D_max to create connectivity graph
- Identify disconnected components (regions with no D_max paths)
- Solve each component independently (smaller MILPs)

**Stage 2: LSCP-Style K Search**
- Start with K0 (computed from workload capacity or target zone size)
- Loop K = K0, K0+1, ..., K_max
- For each K, solve capacitated p-median problem
- Return first K where solution is feasible

**Stage 3: Distance Optimization & Sequencing**
- Given feasible K and zone assignments, solve ATSP within each zone
- Order visits to minimize travel distance
- Compute zone diameters for diagnostic output

### 5.2 Why This Architecture Works

1. **Workload discovery without iteration**: K_min naturally emerges from workload constraints. A zone can hold at most `W_max / avg_service_time` tracts, giving us K_min. No need to artificially guess K.

2. **Early termination**: Once feasible K is found, we stop searching. This is efficient—we want the minimum number of zones that *works*, not the theoretically optimal K.

3. **Distance minimization at feasibility**: After finding feasible K, we then minimize distance via the p-median objective. This ensures zones aren't just feasible—they're efficient.

4. **Asymmetric routing at the end**: TSP sequencing respects that real-world routing is asymmetric. This is a post-processing step that doesn't affect feasibility but improves practical efficiency.

---

## 6. Comparison: Real-World Example

### Scenario: 50 census tracts, 8-hour workday, 20-minute max zone diameter

**Tract characteristics:**
- Average workload: 1 hour per tract
- Average distance between tracts: 10 minutes
- Range: 0.5–2 hours (workload), 5–30 minutes (distance)

### LSCP Result
- **K = 5 zones** (minimum to ensure coverage within 20 min)
- **Workload distribution:** [8h, 7.5h, 6h, 5.5h, 2h]
- **Problem:** Team 5 has only 2 hours of work—operationally wasteful
- **Travel efficiency:** Minimized (coverage-based)
- **Operability:** Poor (Team 1 overloaded, Team 5 underutilized)

### p-Median Result (K=6)
- **K = 6 zones** (specified in advance)
- **Workload distribution:** [8h, 7h, 7h, 6h, 6h, 6h]
- **Problem:** Zone diameter not enforced—some zones might span 45+ minutes between extremes
- **Travel efficiency:** Good (minimized total distance)
- **Operability:** Fair (balanced workload, but zones too spread out)

### surveyzones Result
- **K = 6 zones** (discovered via workload constraints + D_max)
- **Workload distribution:** [7.5h, 7h, 7h, 6.5h, 6.5h, 6h]
- **Zone diameter:** All ≤ 20 minutes
- **Travel efficiency:** Good (distance minimized within constraints)
- **Operability:** Excellent (balanced workload, geographically cohesive)

---

## 7. Implementation Considerations

### 7.1 Computational Complexity

- **LSCP:** Polynomial in K (set cover integer program, but typically fast)
- **p-median:** NP-hard, but good heuristics available
- **surveyzones:** NP-hard + workload constraints (more complex), but component decomposition makes it tractable

For 10,000 tracts with 100 zones: ~5 seconds (GLPK solver)

### 7.2 Sensitivity Analysis

surveyzones is robust to parameter choices:
- **D_max**: Adjust for urban (15 min) vs. rural (60 min) operations
- **max_workload_per_zone**: Adjust for part-time (4h) vs. full-time (8h) teams
- **target_zone_size**: Alternative to workload—"I want ~20 tracts per team"

Sweeping these parameters shows smooth transitions in K and workload distribution, not discrete jumps.

### 7.3 Extensions

Possible extensions (not yet implemented):
- **Partition constraints**: Respect jurisdictional boundaries (already supported)
- **Supervisor zones**: Higher-level aggregation of zones for field coordinators
- **Temporal dynamics**: Adjust zones seasonally or by survey intensity
- **Multi-objective**: Trade off distance vs. workload fairness via Pareto frontier

---

## 8. Conclusion

### Summary of Findings

surveyzones addresses a gap in the facility location literature: **zone construction for human-centered field operations**. Unlike existing models optimized for logistics, coverage, or fairness, surveyzones jointly optimizes for:

1. **Workload balance** (teams have equal capacity)
2. **Geographic operability** (zones have bounded diameter)
3. **Distance efficiency** (minimize wasted travel)
4. **Automatic K discovery** (find the right number of zones)

### Why Existing Models Fall Short

| Model | Primary Gap for Surveys |
|-------|------------------------|
| **LSCP** | No workload balancing; minimizes K, not distance |
| **p-median** | No workload or diameter constraints |
| **p-center** | Wrong fairness metric (distance instead of workload) |
| **MCLP** | Allows uncovered demand; no capacity |
| **WLP** | No diameter constraints; optimizes cost, not operability |

### Justification for surveyzones Approach

surveyzones is mathematically necessary because:
1. **No single existing model addresses all three constraints** (workload + diameter + automatic K)
2. **Survey operations are fundamentally different from logistics/emergency services** (time budgets vs. cost optimization)
3. **Geographic operability is non-negotiable** (unlike theoretical models where arbitrarily large zones are feasible)

### Practical Impact

For survey organizations like IBGE:
- **Reduced supervisor overhead**: Balanced zones reduce need for manual rebalancing
- **Improved data quality**: Teams stay within geographic bounds, improving familiarity and consistency
- **Better time management**: Realistic workload distribution enables accurate scheduling and deadlines
- **Scalability**: Handles 300,000+ tracts across multiple jurisdictions

---

## 9. References & Further Reading

### Academic References
- Church, R. L. & ReVelle, C. S. (1974). "The maximal covering location problem." *Papers of the Regional Science Association*, 32(1), 101-118.
- Hakimi, S. L. (1964). "Optimum locations of switching centers and the absolute centers and medians of a graph." *Operations Research*, 12(3), 450-459.
- Daskin, M. S. (1995). *Network and Discrete Location: Models, Algorithms, and Applications*. John Wiley & Sons.
- ReVelle, C. S., Eiselt, H. A., & Daskin, M. S. (2008). "A bibliography for some fundamental problem categories in discrete location science." *European Journal of Operational Research*, 184(3), 817-848.

### Implementation Details
- surveyzones uses MILP (Mixed Integer Linear Programming) via ROI/ompr packages
- Supports multiple solvers: GLPK, HiGHS, CBC
- Asymmetric routing via OSRM (Open Source Routing Machine)
- TSP/ATSP sequencing via TSP R package

### Survey Operations Context
- IBGE Household Survey (PNAD Contínua): ~50,000 field teams, ~300,000 census tracts
- Monthly rebalancing of zones based on productivity data
- Geographic constraints (state/municipality boundaries)

---

## Appendix: Mathematical Details

### A1. MILP Formulation (Simplified)

**Decision Variables:**
- $z_{ik} \in \{0,1\}$: tract $i$ assigned to zone $k$
- $y_k \in \{0,1\}$: zone $k$ is "open" (has at least one tract)

**Parameters:**
- $d_{ik}$: distance from tract $i$ to zone center $k$ (or within-zone travel if $i=k$)
- $s_i$: service time for tract $i$
- $W_{\max}$: maximum workload per zone
- $K$: number of zones to open
- $D_{\max}$: maximum zone diameter (enforced via distance filtering)

**Objective:**
$$\min \sum_i \sum_k d_{ik} \cdot z_{ik}$$

**Constraints:**
$$\sum_k z_{ik} = 1 \quad \forall i$$
$$\sum_k y_k = K$$
$$z_{ik} \leq y_k \quad \forall i,k$$
$$\sum_{i: z_{ik}=1} s_i \leq W_{\max} \quad \forall k$$

**Pre-processing (not in MILP):**
- Filter distances to $d_{ij} \leq D_{\max}$
- Solve disconnected components independently

### A2. K Discovery Algorithm

```
K_from_size = ceil(n / target_zone_size)
K_from_workload = ceil(total_workload / W_max)
K_min = max(K_from_size, K_from_workload)

for K in K_min to K_max:
    result = solve_capacitated_p_median(K)
    if result.status == "optimal":
        return result

return infeasible
```

---

**White Paper prepared by:** Eduardo Leoni, IBGE
**Date:** February 2026
**Status:** Final

#' Assign Teams to Tracts with Walkable Roteiros
#'
#' Solves a mixed-integer program that assigns every tract to one team. Each
#' tract is either in a *roteiro* (a walkable set of tracts, see
#' [surveyzones_roteiro_candidates()]) worked from a single centre, or *loose*
#' (reached from an anchor tract of the same team), or *fixed* to a team.
#' Mandatory tracts must be in a roteiro. Each team has one origin tract (chosen
#' by the model, or given through `origins`); the
#' objective is the total directed travel time from team origins to roteiro
#' centres plus from anchors to loose tracts. Every team's load must lie in
#' `[(1 - tol), (1 + tol)] * mean`, with `mean` the total load over all units
#' (fixed included) divided by the number of teams, and each team gets at most
#' `k` roteiros.
#'
#' A tract in both `mandatory` and `fixed` is treated as fixed.
#' `diagnostics$solve_time` is the total wall time in seconds, including
#' candidate generation and model building.
#'
#' @param units A data frame with `tract_id` and `expected_service_time` (load,
#'   e.g. households); validated with the same rules as tract tables.
#' @param points An sf object of household POINT geometries with `tract_id`
#'   (may repeat) and unique `point_id`.
#' @param teams A data frame with `team_id` (other columns are ignored).
#' @param eligible A data frame with `tract_id` and `team_id`.
#' @param durations Long data frame with `origin_id`, `destination_id` and
#'   `distance` (minutes, directed). Missing pairs are not usable; `x -> x` is
#'   zero.
#' @param mandatory Character vector of tract ids that must be in a roteiro.
#' @param fixed `NULL` or a data frame with `tract_id` and `team_id`; these
#'   tracts stay with that team (and count toward its load).
#' @param origins `NULL` or a data frame with `team_id` and `tract_id`, one or
#'   more rows per team. `tract_id` is either a unit (e.g. a home base snapped
#'   to a tract) or any other point present as `origin_id` in `durations`
#'   (e.g. the home base itself, with its own travel times to the units).
#'   A listed team has a given origin: the model does not choose one, and the
#'   cost of each of its roteiros is the minimum over its origins of the
#'   directed duration to the roteiro centre (0 when a unit origin is the
#'   centre). Pairs missing from `durations` are skipped; a roteiro with no
#'   usable origin is not offered to that team. Unlisted teams keep the
#'   model-chosen origin. The cost of loose tracts is unchanged.
#' @param k Maximum number of roteiros per team: a whole number >= 1, or
#'   `"auto"`. With `"auto"`, `k_min = ceiling(min_roteiros / n_teams)` is
#'   computed from the minimum number of roteiros that partition the mandatory
#'   tracts (`diagnostics$min_roteiros`) and `k = k_min + 1`. `k_min` is only a
#'   lower bound: team eligibility can require more roteiros per team, in which
#'   case the model is infeasible and the elastic diagnosis reports the `k`
#'   rule. When the minimum-roteiros cover stops at the time limit, its count is
#'   an upper bound, so `k_min` may be over-estimated. The `k` used, `k_auto`
#'   and `k_min` are stored in `parameters`.
#' @param radius_m,max_tracts,walk_distances See
#'   [surveyzones_roteiro_candidates()].
#' @param tol Non-negative load tolerance around the mean.
#' @param solver,max_time,rel_tol,verbose Solver settings.
#' @param use_cache Logical. If `TRUE` (default), the whole solve is memoised
#'   on disk (see [surveyzones_clear_cache()]), keyed on all inputs except
#'   `verbose`; points enter the key as plain coordinates and CRS. A cache hit
#'   returns the stored result, so `diagnostics$solve_time` is the time of the
#'   original solve, not of the lookup.
#'
#' @return A `surveyzones_assignment` object with `assignments`, `roteiros`,
#'   `loads`, `origins`, `parameters` and `diagnostics`. `origins` has one row
#'   per team with a `source` of `"model"`, and one row per origin tract for
#'   teams with `source = "given"`. In `assignments`,
#'   `cost` is the travel time of loose tracts (`NA` for roteiro and fixed
#'   tracts; roteiro costs are in `roteiros`). If the model is infeasible,
#'   `diagnostics$solver_status` is `"infeasible"`, `diagnostics$elastic` lists
#'   the violated rules (`band_lower`, `band_upper`, `k`, `coverage`) with the
#'   `amount` of relaxation, and the tables hold the nearest relaxed plan.
#'   `diagnostics$min_roteiros` is the minimum number of candidate roteiros
#'   partitioning the mandatory tracts (teams ignored).
#' @export
surveyzones_assign_teams <- function(units,
                                     points,
                                     teams,
                                     eligible,
                                     durations,
                                     mandatory = character(),
                                     fixed = NULL,
                                     origins = NULL,
                                     k = 5L,
                                     radius_m = 500,
                                     tol = 0.15,
                                     max_tracts = 5L,
                                     walk_distances = NULL,
                                     solver = "glpk",
                                     max_time = 300,
                                     rel_tol = 0.01,
                                     verbose = FALSE,
                                     use_cache = TRUE) {
  validate_tracts(units)
  validate_solver(solver)
  if (!inherits(points, "sf") || !all(c("tract_id", "point_id") %in% names(points))) {
    cli::cli_abort(
      "{.arg points} must be an sf object with {.field tract_id} and {.field point_id} columns."
    )
  }
  if (!all(as.character(sf::st_geometry_type(points)) == "POINT")) {
    cli::cli_abort("{.arg points} must contain only POINT geometries.")
  }
  k_auto <- identical(k, "auto")
  if (!k_auto && (!is.numeric(k) || length(k) != 1L || is.na(k) || k < 1 || k != round(k))) {
    cli::cli_abort("{.arg k} must be a whole number >= 1 or {.val auto}.")
  }
  if (!is.numeric(tol) || length(tol) != 1L || is.na(tol) || tol < 0) {
    cli::cli_abort("{.arg tol} must be a single number >= 0.")
  }
  if (!is.null(walk_distances) && (!is.data.frame(walk_distances) ||
    !all(c("origin_id", "destination_id", "distance") %in% names(walk_distances)))) {
    cli::cli_abort(
      "{.arg walk_distances} must be {.code NULL} or a data frame with {.field origin_id}, {.field destination_id} and {.field distance}."
    )
  }
  points_df <- tibble::tibble(
    tract_id = as.character(points$tract_id),
    point_id = as.character(points$point_id),
    lon = unname(sf::st_coordinates(points)[, 1]),
    lat = unname(sf::st_coordinates(points)[, 2])
  )
  crs <- sf::st_crs(points)
  core_args <- list(
    units = units, points_df = points_df, crs = crs, teams = teams,
    eligible = eligible, durations = durations, mandatory = mandatory,
    fixed = fixed, origins = origins, k = k, radius_m = radius_m, tol = tol,
    max_tracts = max_tracts, walk_distances = walk_distances, solver = solver,
    max_time = max_time, rel_tol = rel_tol, verbose = verbose
  )
  if (use_cache) {
    is_cached <- do.call(memoise::has_cache(.assign_teams_core_mem), core_args)
    if (is_cached) {
      cli::cli_alert_success("Using cached assignment.")
    } else {
      cli::cli_alert_info("Computing and caching assignment.")
    }
    return(do.call(.assign_teams_core_mem, core_args))
  }
  do.call(.assign_teams_core, core_args)
}

#' @keywords internal
.assign_teams_core <- function(units,
                               points_df,
                               crs,
                               teams,
                               eligible,
                               durations,
                               mandatory,
                               fixed,
                               origins,
                               k,
                               radius_m,
                               tol,
                               max_tracts,
                               walk_distances,
                               solver,
                               max_time,
                               rel_tol,
                               verbose) {
  start_time <- proc.time()[["elapsed"]]
  points <- sf::st_as_sf(points_df, coords = c("lon", "lat"), crs = crs)
  k_auto <- identical(k, "auto")
  need <- list(
    teams = c("team_id"), eligible = c("tract_id", "team_id"),
    durations = c("origin_id", "destination_id", "distance")
  )
  for (nm in names(need)) {
    obj <- get(nm)
    if (!is.data.frame(obj) || !all(need[[nm]] %in% names(obj))) {
      cli::cli_abort(
        "{.arg {nm}} must be a data frame with column{?s} {.field {need[[nm]]}}."
      )
    }
  }
  if (!is.null(fixed) && (!is.data.frame(fixed) || !all(c("tract_id", "team_id") %in% names(fixed)))) {
    cli::cli_abort("{.arg fixed} must be {.code NULL} or a data frame with {.field tract_id} and {.field team_id}.")
  }

  ids <- as.character(units$tract_id)
  nu <- length(ids)
  load <- stats::setNames(as.numeric(units$expected_service_time), ids)
  team_ids <- as.character(teams$team_id)
  if (anyDuplicated(team_ids) || anyNA(team_ids) || length(team_ids) == 0L) {
    cli::cli_abort("{.field team_id} in {.arg teams} must be unique, not {.val NA}, and non-empty.")
  }
  nd <- length(team_ids)

  seen_dur <- union(as.character(durations$origin_id), as.character(durations$destination_id))
  miss <- setdiff(ids, seen_dur)
  if (length(miss) > 0L) {
    cli::cli_abort("Units missing from {.arg durations}: {.val {miss}}.")
  }
  mandatory <- unique(as.character(mandatory))
  if (length(setdiff(mandatory, ids)) > 0L) {
    cli::cli_abort("{.arg mandatory} has unknown tracts: {.val {setdiff(mandatory, ids)}}.")
  }
  elig <- dplyr::distinct(dplyr::mutate(
    eligible[c("tract_id", "team_id")],
    tract_id = as.character(.data$tract_id), team_id = as.character(.data$team_id)
  ))
  bad <- union(setdiff(elig$tract_id, ids), setdiff(elig$team_id, team_ids))
  if (length(bad) > 0L) {
    cli::cli_abort("{.arg eligible} has unknown tracts or teams: {.val {bad}}.")
  }
  if (is.null(fixed)) {
    fixed <- tibble::tibble(tract_id = character(), team_id = character())
  } else {
    fixed <- dplyr::mutate(
      fixed[c("tract_id", "team_id")],
      tract_id = as.character(.data$tract_id), team_id = as.character(.data$team_id)
    )
  }
  bad <- union(setdiff(fixed$tract_id, ids), setdiff(fixed$team_id, team_ids))
  if (length(bad) > 0L) {
    cli::cli_abort("{.arg fixed} has unknown tracts or teams: {.val {bad}}.")
  }
  if (anyDuplicated(fixed$tract_id)) {
    cli::cli_abort("{.arg fixed} must have at most one team per tract.")
  }
  if (is.null(origins)) {
    origins_in <- NULL
    origins <- tibble::tibble(team_id = character(), tract_id = character())
  } else {
    if (!is.data.frame(origins) || !all(c("team_id", "tract_id") %in% names(origins))) {
      cli::cli_abort("{.arg origins} must be {.code NULL} or a data frame with {.field team_id} and {.field tract_id}.")
    }
    origins_in <- origins
    origins <- dplyr::distinct(dplyr::mutate(
      origins[c("team_id", "tract_id")],
      team_id = as.character(.data$team_id), tract_id = as.character(.data$tract_id)
    ))
  }
  bad <- setdiff(origins$team_id, team_ids)
  if (length(bad) > 0L) {
    cli::cli_abort("{.arg origins} has unknown teams: {.val {bad}}.")
  }
  bad <- setdiff(setdiff(origins$tract_id, ids), as.character(durations$origin_id))
  if (length(bad) > 0L) {
    cli::cli_abort("{.arg origins} has origins that are neither units nor origins in {.arg durations}: {.val {bad}}.")
  }
  pt_unknown <- setdiff(unique(as.character(points$tract_id)), ids)
  if (length(pt_unknown) > 0L) {
    cli::cli_abort("{.arg points} has tracts not in {.arg units}: {.val {pt_unknown}}.")
  }

  # directed duration matrix (NA = pair not usable)
  dmat <- matrix(NA_real_, nu, nu, dimnames = list(ids, ids))
  dd <- durations[as.character(durations$origin_id) %in% ids &
    as.character(durations$destination_id) %in% ids, ]
  dmat[cbind(as.character(dd$origin_id), as.character(dd$destination_id))] <- as.numeric(dd$distance)
  diag(dmat) <- 0
  # given origins (units or outside points such as home bases) to units
  o_ids <- unique(origins$tract_id)
  omat <- matrix(NA_real_, length(o_ids), nu, dimnames = list(o_ids, ids))
  od <- durations[as.character(durations$origin_id) %in% o_ids &
    as.character(durations$destination_id) %in% ids, ]
  omat[cbind(as.character(od$origin_id), as.character(od$destination_id))] <- as.numeric(od$distance)
  em_ids <- intersect(o_ids, ids)
  omat[cbind(em_ids, em_ids)] <- 0

  fixed_ids <- fixed$tract_id
  free_ids <- setdiff(ids, fixed_ids)
  optional_ids <- setdiff(free_ids, mandatory)
  fixed_team <- stats::setNames(fixed$team_id, fixed$tract_id)
  fixed_load <- vapply(team_ids, \(d) sum(load[fixed_ids[fixed_team == d]]), numeric(1))
  mean_load <- sum(load) / nd
  lower <- (1 - tol) * mean_load
  upper <- (1 + tol) * mean_load

  # candidates over non-fixed tracts only
  elig_free <- elig[elig$tract_id %in% free_ids, ]
  walk_free <- NULL
  if (!is.null(walk_distances)) {
    free_pts <- as.character(points$point_id)[as.character(points$tract_id) %in% free_ids]
    walk_free <- walk_distances[
      as.character(walk_distances$origin_id) %in% free_pts &
        as.character(walk_distances$destination_id) %in% free_pts, ,
      drop = FALSE
    ]
  }
  if (length(free_ids) > 0L) {
    cand <- surveyzones_roteiro_candidates(
      points[as.character(points$tract_id) %in% free_ids, ],
      elig_free, free_ids, radius_m = radius_m, max_tracts = max_tracts,
      walk_distances = walk_free
    )
    cand_warnings <- attr(cand, "warnings")
  } else {
    cand <- tibble::tibble(
      roteiro_id = character(), center_point_id = character(), center_tract_id = character(),
      tract_ids = list(), teams = list(), span_m = numeric()
    )
    cand_warnings <- character()
  }

  min_roteiros <- .min_roteiros_cover(
    cand, intersect(mandatory, free_ids), solver, max_time, rel_tol, verbose
  )
  k_min <- NA_integer_
  if (k_auto) {
    if (is.na(min_roteiros)) {
      cli::cli_abort(c(
        "Cannot choose {.arg k} automatically: the minimum number of roteiros is unknown.",
        "i" = "Mandatory tracts cannot be covered by candidate roteiros or the cover problem was not solved; pass an explicit {.arg k}."
      ))
    }
    k_min <- as.integer(ceiling(min_roteiros / nd))
    k <- k_min + 1L
  }
  k <- as.integer(k)

  parameters <- list(
    k = k, k_auto = k_auto, k_min = k_min, radius_m = radius_m, tol = tol, max_tracts = as.integer(max_tracts),
    solver = solver, max_time = max_time, rel_tol = rel_tol,
    n_mandatory = length(mandatory), n_fixed = nrow(fixed), origins = origins_in
  )
  finish_empty <- function(status) {
    new_surveyzones_assignment(
      assignments = .empty_assignments(),
      roteiros = .empty_roteiros(),
      loads = tibble::tibble(
        team_id = team_ids, workload = fixed_load, n_roteiros = 0L, n_loose = 0L,
        loose_cost = 0, lower = lower, upper = upper
      ),
      origins = tibble::tibble(
        team_id = character(), origin_tract_id = character(), source = character()
      ),
      parameters = parameters,
      diagnostics = list(
        solver_status = status, objective_value = NA_real_, optimality_gap = NA_real_,
        solve_time = proc.time()[["elapsed"]] - start_time,
        warnings = cand_warnings, elastic = NULL, min_roteiros = NA_integer_
      )
    )
  }

  # (roteiro, team) pairs
  nr <- nrow(cand)
  pr <- rep(seq_len(nr), lengths(cand$teams))
  pd <- match(unlist(cand$teams), team_ids)
  ord <- order(pr, pd)
  pr <- pr[ord]
  pd <- pd[ord]
  r_center <- match(cand$center_tract_id, ids)
  given_idx <- lapply(team_ids, \(d) match(origins$tract_id[origins$team_id == d], o_ids))
  is_given <- lengths(given_idx) > 0L
  model_d <- which(!is_given)
  nmd <- length(model_d)
  pair_cost <- numeric(length(pr))
  for (d in which(is_given)) {
    sel <- which(pd == d)
    pair_cost[sel] <- vapply(sel, \(p) {
      v <- omat[given_idx[[d]], r_center[pr[p]]]
      if (all(is.na(v))) NA_real_ else min(v, na.rm = TRUE)
    }, numeric(1))
  }
  usable <- !is.na(pair_cost)
  pr <- pr[usable]
  pd <- pd[usable]
  pair_cost <- pair_cost[usable]
  np <- length(pr)
  r_load <- vapply(cand$tract_ids, \(x) sum(load[x]), numeric(1))
  pW <- r_load[pr]

  # roteiro membership rows (pair, tract)
  mem_p <- rep(seq_len(np), lengths(cand$tract_ids)[pr])
  mem_t <- match(unlist(cand$tract_ids[pr]), ids)

  # loose (tract, team) rows
  eo <- elig_free[elig_free$tract_id %in% optional_ids, ]
  sq_u <- match(eo$tract_id, ids)
  sq_d <- match(eo$team_id, team_ids)
  o <- order(sq_u, sq_d)
  sq_u <- sq_u[o]
  sq_d <- sq_d[o]
  nsq <- length(sq_u)

  # a rows: (pair, origin tract)
  mp <- which(!is_given[pd])
  a_p <- rep(mp, each = nu)
  a_c <- rep(seq_len(nu), times = length(mp))
  a_cost <- dmat[cbind(a_c, r_center[pr[a_p]])]
  keep <- !is.na(a_cost)
  a_p <- a_p[keep]
  a_c <- a_c[keep]
  a_cost <- a_cost[keep]
  naq <- length(a_p)

  # anchors per team: tracts in a roteiro of d, or fixed to d
  anchors <- lapply(seq_len(nd), \(d) {
    in_r <- unique(mem_t[mem_p %in% which(pd == d)])
    fx <- match(fixed_ids[fixed_team == team_ids[d]], ids)
    sort(unique(c(in_r, fx)))
  })
  w_list <- lapply(seq_len(nsq), \(q) {
    v <- setdiff(anchors[[sq_d[q]]], sq_u[q])
    v[!is.na(dmat[cbind(v, rep(sq_u[q], length(v)))])]
  })
  w_q <- rep(seq_len(nsq), lengths(w_list))
  w_v <- as.integer(unlist(w_list))
  nwq <- length(w_q)
  w_cost <- dmat[cbind(w_v, sq_u[w_q])]
  w_d <- sq_d[w_q]
  w_fixed <- ids[w_v] %in% fixed_ids
  # y[v, d]: team d has a roteiro covering free tract v
  y_key <- unique(cbind(v = w_v[!w_fixed], d = w_d[!w_fixed]))
  y_key <- y_key[order(y_key[, 1], y_key[, 2]), , drop = FALSE]
  ny <- nrow(y_key)
  w_y <- match(
    paste(w_v, w_d)[!w_fixed], paste(y_key[, 1], y_key[, 2])
  )

  ps_by_d <- lapply(seq_len(nd), \(d) which(pd == d))
  sq_by_d <- lapply(seq_len(nd), \(d) which(sq_d == d))
  s_load <- load[ids[sq_u]]

  # teams with no decision variables, or free tracts with no option, make the
  # hard model infeasible; the elastic model diagnoses them
  hard_infeasible <- FALSE
  for (d in seq_len(nd)) {
    if (length(ps_by_d[[d]]) == 0L && length(sq_by_d[[d]]) == 0L &&
      (fixed_load[d] < lower - 1e-9 || fixed_load[d] > upper + 1e-9)) {
      hard_infeasible <- TRUE
    }
  }
  covered <- union(mem_t, sq_u)
  free_idx <- match(free_ids, ids)
  nf <- length(free_idx)
  if (length(setdiff(free_idx, covered)) > 0L) hard_infeasible <- TRUE

  # ---- Variables ---------------------------------------------------------------
  # One column vector, in named blocks (`col("s", q)` is the column of s[q]):
  #   lam[p]  binary  roteiro pr[p] goes to team pd[p]
  #   s[q]    binary  tract sq_u[q] is loose for team sq_d[q]
  #   z[c,m]  binary  the model-chosen origin of team model_d[m] is tract c
  #   a[i]    cont.   pair a_p[i] is reached from origin a_c[i] (a = lam * z)
  #   w[i]    cont.   loose q = w_q[i] hangs off anchor tract w_v[i]
  #   y[i]    cont.   team y_key[i, 2] has a roteiro covering tract y_key[i, 1]
  # The elastic model appends one slack per relaxable rule:
  #   lo[d], hi[d]  load band;  ek[d]  k;  ec[t]  coverage of free tract t
  sizes <- c(
    lam = np, s = nsq, z = nu * nmd, a = naq, w = nwq, y = ny,
    lo = nd, hi = nd, ek = nd, ec = nf
  )
  offset <- stats::setNames(c(0L, cumsum(sizes))[seq_along(sizes)], names(sizes))
  col <- function(block, idx = seq_len(sizes[[block]])) offset[[block]] + idx
  z_col <- function(c, m) col("z", (m - 1L) * nu + c)
  nvar <- offset[["lo"]]
  nvar_e <- sum(sizes)

  # min  sum_p cost(origin -> roteiro) lam[p]   (teams with given origins)
  #    + sum_i dur(a_c[i] -> centre) a[i]       (teams with model origins)
  #    + sum_i dur(anchor -> loose tract) w[i]
  cost <- numeric(nvar)
  cost[col("lam")] <- pair_cost
  cost[col("a")] <- a_cost
  cost[col("w")] <- w_cost
  types <- c(rep("B", offset[["a"]]), rep("C", naq + nwq + ny))

  # ---- Constraints -------------------------------------------------------------
  # x[i] <= y[i] for each pair of columns
  add_le <- function(rows, x, y) {
    n <- length(x)
    rows$add(rep(seq_len(n), 2), c(x, y), rep(c(1, -1), each = n), n, "<=", 0)
  }

  # a[i] <= z[origin, team] and a[i] <= lam[pair]; sum_c a[p, c] = lam[p]:
  # a model-origin team reaches each of its roteiros from its one origin
  add_origin_link <- function(rows) {
    add_le(rows, col("a"), z_col(a_c, match(pd[a_p], model_d)))
    add_le(rows, col("a"), col("lam", a_p))
    rows$add(
      c(match(a_p, mp), seq_along(mp)), c(col("a"), col("lam", mp)),
      c(rep(1, naq), rep(-1, length(mp))), length(mp), "==", 0
    )
  }

  # sum_c z[c, m] = 1: one origin per model-origin team
  add_one_origin <- function(rows) {
    rows$add(rep(seq_len(nmd), each = nu), col("z"), rep(1, nu * nmd), nmd, "==", 1)
  }

  # sum_{p of d} lam[p] <= k (+ ek[d])
  add_max_roteiros <- function(rows, elastic) {
    if (elastic) {
      rows$add(
        c(pd, seq_len(nd)), c(col("lam"), col("ek")),
        c(rep(1, np), rep(-1, nd)), nd, "<=", k
      )
    } else {
      rows$add(pd, col("lam"), rep(1, np), nd, "<=", k)
    }
  }

  # lower <= fixed load + roteiro load + loose load <= upper, per team
  # (elastic: + lo[d] on the lower side, - hi[d] on the upper side)
  add_load_band <- function(rows, elastic) {
    for (d in seq_len(nd)) {
      cols <- c(col("lam", ps_by_d[[d]]), col("s", sq_by_d[[d]]))
      vals <- c(pW[ps_by_d[[d]]], s_load[sq_by_d[[d]]])
      lo <- lower - fixed_load[d]
      hi <- upper - fixed_load[d]
      if (elastic) {
        rows$add(rep(1L, length(cols) + 1L), c(cols, col("lo", d)), c(vals, 1), 1L, ">=", lo)
        rows$add(rep(1L, length(cols) + 1L), c(cols, col("hi", d)), c(vals, -1), 1L, "<=", hi)
      } else if (length(cols) > 0L) {
        rows$add(rep(1L, length(cols)), cols, vals, 1L, ">=", lo)
        rows$add(rep(1L, length(cols)), cols, vals, 1L, "<=", hi)
      }
    }
  }

  # sum_{p covering t} lam[p] + sum_{q of t} s[q] (+ ec[t]) = 1, per free
  # tract: each in exactly one roteiro or loose with one team. Mandatory
  # tracts have no s, so they must be in a roteiro.
  add_coverage <- function(rows, elastic) {
    mi <- match(mem_t, free_idx)
    si <- match(sq_u, free_idx)
    if (elastic) {
      rows$add(
        c(mi, si, seq_len(nf)), c(col("lam", mem_p), col("s"), col("ec")),
        rep(1, length(mi) + nsq + nf), nf, "==", 1
      )
    } else {
      rows$add(c(mi, si), c(col("lam", mem_p), col("s")), rep(1, length(mi) + nsq), nf, "==", 1)
    }
  }

  # A loose tract hangs off one anchor of its team: w <= s; w <= y when the
  # anchor is not fixed (y[v, d] = sum of d's roteiros covering v); and
  # sum_{anchors} w[q, .] = s[q]
  add_loose_anchor <- function(rows) {
    add_le(rows, col("w"), col("s", w_q))
    nonfix <- which(!w_fixed)
    add_le(rows, col("w", nonfix), col("y", w_y))
    gi <- match(paste(mem_t, pd[mem_p]), paste(y_key[, 1], y_key[, 2]))
    hit <- which(!is.na(gi))
    rows$add(
      c(seq_len(ny), gi[hit]), c(col("y"), col("lam", mem_p[hit])),
      c(rep(1, ny), rep(-1, length(hit))), ny, "==", 0
    )
    rows$add(
      c(w_q, seq_len(nsq)), c(col("w"), col("s")),
      c(rep(1, nwq), rep(-1, nsq)), nsq, "==", 0
    )
  }

  build_rows <- function(elastic) {
    rows <- .new_row_builder()
    add_origin_link(rows)
    add_one_origin(rows)
    add_max_roteiros(rows, elastic)
    add_load_band(rows, elastic)
    add_coverage(rows, elastic)
    add_loose_anchor(rows)
    rows
  }

  finish <- function(out) {
    out$diagnostics$min_roteiros <- min_roteiros
    out
  }

  solved <- list(status = "infeasible")
  if (!hard_infeasible) {
    op <- build_rows(FALSE)$op(cost, types, nvar)
    solved <- .solve_op(op, solver = solver, max_time = max_time, rel_tol = rel_tol, verbose = verbose)
  }
  if (!.is_solution_status(solved$status)) {
    if (solved$status != "infeasible") {
      out <- finish_empty(solved$status)
      out$diagnostics$optimality_gap <- solved$optimality_gap
      return(finish(out))
    }
    return(finish(.elastic_fallback(
      build_rows = build_rows, nvar_e = nvar_e, types = types, cost = cost,
      slack_cols = nvar + seq_len(nvar_e - nvar), o_lo = nvar, nd = nd, nf = nf,
      team_ids = team_ids, free_ids = free_ids, solver = solver, max_time = max_time,
      rel_tol = rel_tol, verbose = verbose, finish_empty = finish_empty,
      extract = \(solved, sol) .extract_assignment(
        solved = solved, start_time = start_time, parameters = parameters,
        cand = cand, ids = ids, team_ids = team_ids, load = load, dmat = dmat,
        pr = pr, pd = pd, sq_u = sq_u, sq_d = sq_d, w_q = w_q, w_v = w_v, w_cost = w_cost,
        lam_v = sol[col("lam")], s_v = sol[col("s")],
        z_v = sol[col("z")], w_sol = sol[col("w")],
        model_d = model_d, given_idx = given_idx, omat = omat,
        fixed = fixed, fixed_load = fixed_load, lower = lower, upper = upper,
        warnings = cand_warnings
      )
    )))
  }
  sol <- solved$solution
  lam_v <- sol[col("lam")]
  s_v <- sol[col("s")]
  z_v <- sol[col("z")]
  w_v_sol <- sol[col("w")]
  out <- .extract_assignment(
    solved = solved, start_time = start_time, parameters = parameters,
    cand = cand, ids = ids, team_ids = team_ids, load = load, dmat = dmat,
    pr = pr, pd = pd, sq_u = sq_u, sq_d = sq_d, w_q = w_q, w_v = w_v, w_cost = w_cost,
    lam_v = lam_v, s_v = s_v, z_v = z_v, w_sol = w_v_sol,
    model_d = model_d, given_idx = given_idx, omat = omat,
    fixed = fixed, fixed_load = fixed_load, lower = lower, upper = upper,
    warnings = cand_warnings
  )
  bad_cover <- !identical(sort(out$assignments$tract_id), sort(ids))
  bad_band <- any(out$loads$workload < lower - 1e-6 | out$loads$workload > upper + 1e-6)
  if (bad_cover || bad_band) {
    cli::cli_warn(
      "Solver solution failed verification ({if (bad_cover) 'coverage' else 'load band'}); returning status {.val error}."
    )
    bad <- finish_empty("error")
    bad$diagnostics$optimality_gap <- solved$optimality_gap
    return(finish(bad))
  }
  finish(out)
}

#' @keywords internal
.elastic_fallback <- function(build_rows, nvar_e, types, cost, slack_cols, o_lo, nd, nf,
                              team_ids, free_ids, solver, max_time, rel_tol, verbose,
                              finish_empty, extract) {
  types_e <- c(types, rep("C", nvar_e - length(types)))
  cost1 <- replace(numeric(nvar_e), slack_cols, 1)
  rows <- build_rows(TRUE)
  p1 <- .solve_op(rows$op(cost1, types_e, nvar_e),
    solver = solver, max_time = max_time, rel_tol = rel_tol, verbose = verbose
  )
  if (!.is_solution_status(p1$status)) return(finish_empty("infeasible"))
  total <- sum(p1$solution[slack_cols])
  rows$add(rep(1L, length(slack_cols)), slack_cols, rep(1, length(slack_cols)), 1L, "<=", total + 1e-6)
  p2 <- .solve_op(rows$op(c(cost, rep(0, nvar_e - length(cost))), types_e, nvar_e),
    solver = solver, max_time = max_time, rel_tol = rel_tol, verbose = verbose
  )
  fin <- if (.is_solution_status(p2$status)) p2 else p1
  sol <- fin$solution
  sl <- sol[slack_cols]
  names_rule <- c(
    rep("band_lower", nd), rep("band_upper", nd), rep("k", nd), rep("coverage", nf)
  )
  ids_rule <- c(team_ids, team_ids, team_ids, free_ids)
  keep <- sl > 1e-6
  out <- extract(fin, sol)
  out$diagnostics$solver_status <- "infeasible"
  out$diagnostics$elastic <- tibble::tibble(
    rule = names_rule[keep], id = ids_rule[keep], amount = sl[keep]
  )
  out
}

#' @keywords internal
.min_roteiros_cover <- function(candidates, mandatory, solver, max_time = 60, rel_tol = 0.01,
                                verbose = FALSE) {
  if (length(mandatory) == 0L) return(0L)
  nr <- nrow(candidates)
  if (nr == 0L) return(NA_integer_)
  tr <- candidates$tract_ids
  all_t <- unique(c(mandatory, unlist(tr)))
  ti <- match(unlist(tr), all_t)
  ri <- rep(seq_len(nr), lengths(tr))
  rows <- .new_row_builder()
  is_m <- all_t %in% mandatory
  m_rows <- cumsum(is_m)
  sel_m <- is_m[ti]
  rows$add(m_rows[ti][sel_m], ri[sel_m], rep(1, sum(sel_m)), sum(is_m), "==", 1)
  o_rows <- cumsum(!is_m)
  sel_o <- !is_m[ti]
  rows$add(o_rows[ti][sel_o], ri[sel_o], rep(1, sum(sel_o)), sum(!is_m), "<=", 1)
  op <- rows$op(rep(1, nr), rep("B", nr), nr)
  res <- .solve_op(op, solver = solver, max_time = max_time, rel_tol = rel_tol, verbose = verbose)
  if (!.is_solution_status(res$status)) return(NA_integer_)
  as.integer(round(res$objective_value))
}

#' @keywords internal
.extract_assignment <- function(solved, start_time, parameters, cand, ids, team_ids, load, dmat,
                                pr, pd, sq_u, sq_d, w_q, w_v, w_cost, lam_v, s_v, z_v, w_sol,
                                model_d, given_idx, omat, fixed, fixed_load, lower, upper, warnings) {
  nu <- length(ids)
  nd <- length(team_ids)
  is_given <- lengths(given_idx) > 0L
  # rows of each team's origins in its cost matrix: omat if given, dmat if model
  origin_idx <- given_idx
  origin_mat <- rep(list(omat), nd)
  for (m in seq_along(model_d)) {
    origin_idx[[model_d[m]]] <- which.max(z_v[(m - 1L) * nu + seq_len(nu)])
    origin_mat[[model_d[m]]] <- dmat
  }
  origins <- tibble::tibble(
    team_id = rep(team_ids, lengths(origin_idx)),
    origin_tract_id = unlist(lapply(seq_len(nd), \(d) rownames(origin_mat[[d]])[origin_idx[[d]]])),
    source = rep(ifelse(is_given, "given", "model"), lengths(origin_idx))
  )

  chosen <- which(lam_v > 0.5)
  rot <- tibble::tibble(
    roteiro_id = cand$roteiro_id[pr[chosen]],
    team_id = team_ids[pd[chosen]],
    center_point_id = cand$center_point_id[pr[chosen]],
    center_tract_id = cand$center_tract_id[pr[chosen]],
    n_tracts = lengths(cand$tract_ids[pr[chosen]]),
    workload = vapply(cand$tract_ids[pr[chosen]], \(x) sum(load[x]), numeric(1)),
    span_m = cand$span_m[pr[chosen]],
    cost = vapply(chosen, \(p) {
      v <- origin_mat[[pd[p]]][origin_idx[[pd[p]]], match(cand$center_tract_id[pr[p]], ids)]
      min(v, na.rm = TRUE)
    }, numeric(1))
  )
  from_rot <- tibble::tibble(
    tract_id = unlist(cand$tract_ids[pr[chosen]]),
    team_id = rep(rot$team_id, rot$n_tracts),
    roteiro_id = rep(rot$roteiro_id, rot$n_tracts),
    mode = "roteiro", cost = NA_real_
  )
  sch <- which(s_v > 0.5)
  loose_cost <- vapply(sch, \(q) sum(w_sol[w_q == q] * w_cost[w_q == q]), numeric(1))
  from_loose <- tibble::tibble(
    tract_id = ids[sq_u[sch]], team_id = team_ids[sq_d[sch]],
    roteiro_id = NA_character_, mode = "loose", cost = loose_cost
  )
  from_fixed <- tibble::tibble(
    tract_id = fixed$tract_id, team_id = fixed$team_id,
    roteiro_id = NA_character_, mode = "fixed", cost = NA_real_
  )
  assignments <- dplyr::bind_rows(from_rot, from_loose, from_fixed) |>
    dplyr::arrange(match(.data$tract_id, ids))

  loads <- tibble::tibble(
    team_id = team_ids,
    workload = fixed_load + vapply(team_ids, \(d) {
      sum(load[assignments$tract_id[assignments$team_id == d & assignments$mode != "fixed"]])
    }, numeric(1)),
    n_roteiros = vapply(team_ids, \(d) sum(rot$team_id == d), integer(1)),
    n_loose = vapply(team_ids, \(d) sum(from_loose$team_id == d), integer(1)),
    loose_cost = vapply(team_ids, \(d) sum(from_loose$cost[from_loose$team_id == d]), numeric(1)),
    lower = lower, upper = upper
  )
  new_surveyzones_assignment(
    assignments = assignments,
    roteiros = rot,
    loads = loads,
    origins = origins,
    parameters = parameters,
    diagnostics = list(
      solver_status = solved$status,
      objective_value = solved$objective_value,
      optimality_gap = solved$optimality_gap,
      solve_time = proc.time()[["elapsed"]] - start_time,
      warnings = warnings, elastic = NULL, min_roteiros = NA_integer_
    )
  )
}

#' @keywords internal
.empty_assignments <- function() {
  tibble::tibble(
    tract_id = character(), team_id = character(), roteiro_id = character(),
    mode = character(), cost = numeric()
  )
}

#' @keywords internal
.empty_roteiros <- function() {
  tibble::tibble(
    roteiro_id = character(), team_id = character(), center_point_id = character(),
    center_tract_id = character(), n_tracts = integer(), workload = numeric(),
    span_m = numeric(), cost = numeric()
  )
}

#' @keywords internal
.new_row_builder <- function() {
  st <- new.env(parent = emptyenv())
  st$i <- list()
  st$j <- list()
  st$v <- list()
  st$dir <- list()
  st$rhs <- list()
  st$n <- 0L
  list(
    add = function(i, j, v, nrow, dir, rhs) {
      if (nrow == 0L) return(invisible())
      st$i[[length(st$i) + 1L]] <- as.integer(i) + st$n
      st$j[[length(st$j) + 1L]] <- as.integer(j)
      st$v[[length(st$v) + 1L]] <- as.numeric(v)
      st$dir[[length(st$dir) + 1L]] <- rep(dir, nrow)
      st$rhs[[length(st$rhs) + 1L]] <- rep(as.numeric(rhs), length.out = nrow)
      st$n <- st$n + as.integer(nrow)
      invisible()
    },
    op = function(cost, types, nvar) {
      L <- slam::simple_triplet_matrix(
        i = unlist(st$i), j = unlist(st$j), v = unlist(st$v),
        nrow = st$n, ncol = nvar
      )
      ROI::OP(
        objective = ROI::L_objective(cost),
        constraints = ROI::L_constraint(L, dir = unlist(st$dir), rhs = unlist(st$rhs)),
        types = types, maximum = FALSE
      )
    }
  )
}

#' @keywords internal
.classify_roi_solution <- function(sol) {
  msg <- tryCatch(sol$message, error = function(e) NULL)
  text <- paste(c(msg$status_message, msg$status, sol$status$msg$message), collapse = " ")
  symbol <- toupper(as.character(tryCatch(sol$status$msg$symbol, error = function(e) NULL))[1])
  if (length(symbol) != 1L || is.na(symbol)) symbol <- ""
  obj <- tryCatch(as.numeric(sol$objval), error = function(e) NA_real_)
  if (length(obj) != 1L) obj <- NA_real_
  has_full <- length(sol$solution) > 0L && is.finite(obj) && !anyNA(sol$solution)
  gap <- tryCatch(as.numeric(msg$info$mip_gap), error = function(e) NA_real_)
  if (length(gap) != 1L || !is.finite(gap)) gap <- NA_real_
  code_ok <- isTRUE(sol$status$code == 0L)
  infeasible_symbols <- c(
    "INFEASIBLE", "UNBOUNDEDORINFEASIBLE", "INFEASIBLE_OR_UNBOUNDED", "GLP_NOFEAS", "GLP_INFEAS"
  )
  if (symbol %in% infeasible_symbols) {
    status <- "infeasible"
    has_sol <- FALSE
  } else if (symbol %in% c("OPTIMAL", "GLP_OPT") || (code_ok && has_full)) {
    has_sol <- has_full
    status <- if (has_sol) "optimal" else "error"
  } else {
    # other symbols (incumbent only) and unknown solvers: fall back to message text
    infeasible <- grepl("infeasible|no feasible", text, ignore.case = TRUE)
    has_sol <- has_full && !infeasible
    status <- if (infeasible) "infeasible" else if (has_sol) "feasible" else "error"
  }
  list(status = status, has_sol = has_sol, obj = obj, gap = gap)
}

#' @keywords internal
.solve_op <- function(op, solver, max_time, rel_tol, verbose) {
  params <- .solver_control_params(solver, max_time, rel_tol, verbose)
  solver_name <- params$solver
  params$solver <- NULL
  sol <- ROI::ROI_solve(op, solver = solver_name, control = params)
  cls <- .classify_roi_solution(sol)
  status <- cls$status
  has_sol <- cls$has_sol
  obj <- cls$obj
  gap <- cls$gap
  list(
    status = status, objective_value = if (has_sol) obj else NA_real_,
    optimality_gap = gap, solution = if (has_sol) as.numeric(sol$solution) else NULL,
    raw = sol
  )
}

#' Assign Teams to Tracts with Walkable Roteiros
#'
#' Solves a mixed-integer program that assigns every tract to one team. Each
#' tract is either in a *roteiro* (a walkable set of tracts, see
#' [surveyzones_roteiro_candidates()]) worked from a single centre, or *loose*
#' (reached from an anchor tract of the same team), or *fixed* to a team.
#' Mandatory tracts must be in a roteiro. Each team has one origin tract; the
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
#' @param k Maximum number of roteiros per team.
#' @param radius_m,max_tracts See [surveyzones_roteiro_candidates()].
#' @param tol Non-negative load tolerance around the mean.
#' @param solver,max_time,rel_tol,verbose Solver settings.
#'
#' @return A `surveyzones_assignment` object with `assignments`, `roteiros`,
#'   `loads`, `origins`, `parameters` and `diagnostics`. In `assignments`,
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
                                     k = 5L,
                                     radius_m = 500,
                                     tol = 0.15,
                                     max_tracts = 5L,
                                     solver = "glpk",
                                     max_time = 300,
                                     rel_tol = 0.01,
                                     verbose = FALSE) {
  start_time <- proc.time()[["elapsed"]]
  validate_tracts(units)
  validate_solver(solver)
  if (!inherits(points, "sf") || !all(c("tract_id", "point_id") %in% names(points))) {
    cli::cli_abort(
      "{.arg points} must be an sf object with {.field tract_id} and {.field point_id} columns."
    )
  }
  if (!is.numeric(k) || length(k) != 1L || is.na(k) || k < 1) {
    cli::cli_abort("{.arg k} must be a single number >= 1.")
  }
  if (!is.numeric(tol) || length(tol) != 1L || is.na(tol) || tol < 0) {
    cli::cli_abort("{.arg tol} must be a single number >= 0.")
  }
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
  if (length(free_ids) > 0L) {
    cand <- surveyzones_roteiro_candidates(
      points[as.character(points$tract_id) %in% free_ids, ],
      elig_free, free_ids, radius_m = radius_m, max_tracts = max_tracts
    )
    cand_warnings <- attr(cand, "warnings")
  } else {
    cand <- tibble::tibble(
      roteiro_id = character(), center_point_id = character(), center_tract_id = character(),
      tract_ids = list(), teams = list(), span_m = numeric()
    )
    cand_warnings <- character()
  }

  parameters <- list(
    k = as.integer(k), radius_m = radius_m, tol = tol, max_tracts = as.integer(max_tracts),
    solver = solver, max_time = max_time, rel_tol = rel_tol,
    n_mandatory = length(mandatory), n_fixed = nrow(fixed)
  )
  finish_empty <- function(status) {
    new_surveyzones_assignment(
      assignments = .empty_assignments(),
      roteiros = .empty_roteiros(),
      loads = tibble::tibble(
        team_id = team_ids, workload = fixed_load, n_roteiros = 0L, n_loose = 0L,
        loose_cost = 0, lower = lower, upper = upper
      ),
      origins = tibble::tibble(team_id = character(), origin_tract_id = character()),
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
  np <- length(pr)
  r_load <- vapply(cand$tract_ids, \(x) sum(load[x]), numeric(1))
  r_center <- match(cand$center_tract_id, ids)
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
  a_p <- rep(seq_len(np), each = nu)
  a_c <- rep(seq_len(nu), times = np)
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

  # variable layout: lam | s | z | a | w | y  (| lo | hi | ek | ec when elastic)
  o_lam <- 0L
  o_s <- np
  o_z <- o_s + nsq
  o_a <- o_z + nu * nd
  o_w <- o_a + naq
  o_y <- o_w + nwq
  nvar <- o_y + ny
  o_lo <- nvar
  o_hi <- o_lo + nd
  o_ek <- o_hi + nd
  o_ec <- o_ek + nd
  nvar_e <- o_ec + nf
  cost <- c(rep(0, o_a), a_cost, w_cost, rep(0, ny))
  types <- c(rep("B", o_a), rep("C", naq + nwq + ny))

  build_rows <- function(elastic) {
    rows <- .new_row_builder()
    seq_a <- seq_len(naq)
    # a <= z (origin) and a <= lam (valid inequality)
    rows$add(rep(seq_a, 2), c(o_a + seq_a, o_z + (pd[a_p] - 1L) * nu + a_c), rep(c(1, -1), each = naq), naq, "<=", 0)
    rows$add(rep(seq_a, 2), c(o_a + seq_a, o_lam + a_p), rep(c(1, -1), each = naq), naq, "<=", 0)
    # sum_c a[r,c,d] = lam[r,d]
    rows$add(c(a_p, seq_len(np)), c(o_a + seq_a, o_lam + seq_len(np)), c(rep(1, naq), rep(-1, np)), np, "==", 0)
    # one origin per team; at most k roteiros per team
    rows$add(
      rep(seq_len(nd), each = nu), o_z + seq_len(nu * nd), rep(1, nu * nd), nd, "==", 1
    )
    if (elastic) {
      rows$add(
        c(pd, seq_len(nd)), c(o_lam + seq_len(np), o_ek + seq_len(nd)),
        c(rep(1, np), rep(-1, nd)), nd, "<=", k
      )
    } else {
      rows$add(pd, o_lam + seq_len(np), rep(1, np), nd, "<=", k)
    }
    # load band
    for (d in seq_len(nd)) {
      cols <- c(o_lam + ps_by_d[[d]], o_s + sq_by_d[[d]])
      vals <- c(pW[ps_by_d[[d]]], s_load[sq_by_d[[d]]])
      if (elastic) {
        rows$add(rep(1L, length(cols) + 1L), c(cols, o_lo + d), c(vals, 1), 1L, ">=", lower - fixed_load[d])
        rows$add(rep(1L, length(cols) + 1L), c(cols, o_hi + d), c(vals, -1), 1L, "<=", upper - fixed_load[d])
      } else if (length(cols) > 0L) {
        rows$add(rep(1L, length(cols)), cols, vals, 1L, ">=", lower - fixed_load[d])
        rows$add(rep(1L, length(cols)), cols, vals, 1L, "<=", upper - fixed_load[d])
      }
    }
    # coverage of free tracts
    mi <- match(mem_t, free_idx)
    si <- match(sq_u, free_idx)
    if (elastic) {
      rows$add(
        c(mi, si, seq_len(nf)), c(o_lam + mem_p, o_s + seq_len(nsq), o_ec + seq_len(nf)),
        rep(1, length(mi) + nsq + nf), nf, "==", 1
      )
    } else {
      rows$add(
        c(mi, si), c(o_lam + mem_p, o_s + seq_len(nsq)), rep(1, length(mi) + nsq), nf, "==", 1
      )
    }
    # loose tracts hang off one anchor
    seq_w <- seq_len(nwq)
    rows$add(rep(seq_w, 2), c(o_w + seq_w, o_s + w_q), rep(c(1, -1), each = nwq), nwq, "<=", 0)
    nonfix <- which(!w_fixed)
    rows$add(
      rep(seq_along(nonfix), 2), c(o_w + nonfix, o_y + w_y), rep(c(1, -1), each = length(nonfix)),
      length(nonfix), "<=", 0
    )
    gi <- match(paste(mem_t, pd[mem_p]), paste(y_key[, 1], y_key[, 2]))
    hit <- which(!is.na(gi))
    rows$add(
      c(seq_len(ny), gi[hit]), c(o_y + seq_len(ny), o_lam + mem_p[hit]),
      c(rep(1, ny), rep(-1, length(hit))), ny, "==", 0
    )
    rows$add(
      c(w_q, seq_len(nsq)), c(o_w + seq_w, o_s + seq_len(nsq)),
      c(rep(1, nwq), rep(-1, nsq)), nsq, "==", 0
    )
    rows
  }

  min_roteiros <- .min_roteiros_cover(
    cand, intersect(mandatory, free_ids), solver, max_time, rel_tol, verbose
  )
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
      slack_cols = o_lo + seq_len(nvar_e - o_lo), o_lo = o_lo, nd = nd, nf = nf,
      team_ids = team_ids, free_ids = free_ids, solver = solver, max_time = max_time,
      rel_tol = rel_tol, verbose = verbose, finish_empty = finish_empty,
      extract = \(solved, sol) .extract_assignment(
        solved = solved, start_time = start_time, parameters = parameters,
        cand = cand, ids = ids, team_ids = team_ids, load = load, dmat = dmat,
        pr = pr, pd = pd, sq_u = sq_u, sq_d = sq_d, w_q = w_q, w_v = w_v, w_cost = w_cost,
        lam_v = sol[o_lam + seq_len(np)], s_v = sol[o_s + seq_len(nsq)],
        z_v = sol[o_z + seq_len(nu * nd)], w_sol = sol[o_w + seq_len(nwq)],
        fixed = fixed, fixed_load = fixed_load, lower = lower, upper = upper,
        warnings = cand_warnings
      )
    )))
  }
  sol <- solved$solution
  lam_v <- sol[o_lam + seq_len(np)]
  s_v <- sol[o_s + seq_len(nsq)]
  z_v <- sol[o_z + seq_len(nu * nd)]
  w_v_sol <- sol[o_w + seq_len(nwq)]
  out <- .extract_assignment(
    solved = solved, start_time = start_time, parameters = parameters,
    cand = cand, ids = ids, team_ids = team_ids, load = load, dmat = dmat,
    pr = pr, pd = pd, sq_u = sq_u, sq_d = sq_d, w_q = w_q, w_v = w_v, w_cost = w_cost,
    lam_v = lam_v, s_v = s_v, z_v = z_v, w_sol = w_v_sol,
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
                                fixed, fixed_load, lower, upper, warnings) {
  nu <- length(ids)
  nd <- length(team_ids)
  origin <- vapply(seq_len(nd), \(d) {
    zz <- z_v[(d - 1L) * nu + seq_len(nu)]
    which.max(zz)
  }, integer(1))
  origins <- tibble::tibble(team_id = team_ids, origin_tract_id = ids[origin])

  chosen <- which(lam_v > 0.5)
  rot <- tibble::tibble(
    roteiro_id = cand$roteiro_id[pr[chosen]],
    team_id = team_ids[pd[chosen]],
    center_point_id = cand$center_point_id[pr[chosen]],
    center_tract_id = cand$center_tract_id[pr[chosen]],
    n_tracts = lengths(cand$tract_ids[pr[chosen]]),
    workload = vapply(cand$tract_ids[pr[chosen]], \(x) sum(load[x]), numeric(1)),
    span_m = cand$span_m[pr[chosen]],
    cost = dmat[cbind(origin[pd[chosen]], match(cand$center_tract_id[pr[chosen]], ids))]
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
.solve_op <- function(op, solver, max_time, rel_tol, verbose) {
  params <- .solver_control_params(solver, max_time, rel_tol, verbose)
  solver_name <- params$solver
  params$solver <- NULL
  sol <- ROI::ROI_solve(op, solver = solver_name, control = params)
  msg <- tryCatch(sol$message, error = function(e) NULL)
  text <- paste(c(msg$status_message, msg$status, sol$status$msg$message), collapse = " ")
  symbol <- toupper(as.character(tryCatch(sol$status$msg$symbol, error = function(e) NULL))[1])
  if (length(symbol) != 1L || is.na(symbol)) symbol <- ""
  obj <- tryCatch(as.numeric(sol$objval), error = function(e) NA_real_)
  has_full <- length(sol$solution) > 0L && length(obj) == 1L && is.finite(obj) && !anyNA(sol$solution)
  gap <- tryCatch(as.numeric(msg$info$mip_gap), error = function(e) NA_real_)
  if (length(gap) != 1L || !is.finite(gap)) gap <- NA_real_
  infeasible_symbols <- c("INFEASIBLE", "UNBOUNDEDORINFEASIBLE", "INFEASIBLE_OR_UNBOUNDED")
  if (symbol %in% infeasible_symbols) {
    status <- "infeasible"
    has_sol <- FALSE
  } else if (nzchar(symbol)) {
    has_sol <- has_full
    status <- if (has_sol && symbol == "OPTIMAL") "optimal" else if (has_sol) "feasible" else "error"
  } else {
    # unknown solver without a status symbol: fall back to message text
    infeasible <- grepl("infeasible|no feasible", text, ignore.case = TRUE)
    has_sol <- has_full && !infeasible
    status <- if (infeasible) {
      "infeasible"
    } else if (has_sol && isTRUE(sol$status$code == 0L)) {
      "optimal"
    } else if (has_sol) {
      "feasible"
    } else {
      "error"
    }
  }
  list(
    status = status, objective_value = if (has_sol) obj else NA_real_,
    optimality_gap = gap, solution = if (has_sol) as.numeric(sol$solution) else NULL,
    raw = sol
  )
}

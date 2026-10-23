utils::globalVariables(c("lam", "s", "z", "a", "w", "y", "p", "q", "c", "i", "m", "g", "j"))

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
#'   tracts; roteiro costs are in `roteiros`). If the model is infeasible the
#'   tables are empty and `diagnostics$solver_status` is `"infeasible"`.
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
  w_q <- integer()
  w_v <- integer()
  for (q in seq_len(nsq)) {
    v <- setdiff(anchors[[sq_d[q]]], sq_u[q])
    v <- v[!is.na(dmat[cbind(v, sq_u[q])])]
    w_q <- c(w_q, rep(q, length(v)))
    w_v <- c(w_v, v)
  }
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

  # teams with no decision variables: only fixed load matters
  for (d in seq_len(nd)) {
    if (length(ps_by_d[[d]]) == 0L && length(sq_by_d[[d]]) == 0L &&
      (fixed_load[d] < lower - 1e-9 || fixed_load[d] > upper + 1e-9)) {
      return(finish_empty("infeasible"))
    }
  }
  # coverage needs at least one option per free tract
  covered <- union(mem_t, sq_u)
  if (length(setdiff(match(free_ids, ids), covered)) > 0L) {
    return(finish_empty("infeasible"))
  }

  n1 <- \(n) n + 1L
  model <- ompr::MILPModel() |>
    ompr::add_variable(lam[p], p = seq_len(n1(np)), type = "binary") |>
    ompr::add_variable(s[q], q = seq_len(n1(nsq)), type = "binary") |>
    ompr::add_variable(z[c], c = seq_len(nu * nd), type = "binary") |>
    ompr::add_variable(a[i], i = seq_len(n1(naq)), type = "continuous", lb = 0) |>
    ompr::add_variable(w[m], m = seq_len(n1(nwq)), type = "continuous", lb = 0) |>
    ompr::add_variable(y[g], g = seq_len(n1(ny)), type = "continuous", lb = 0) |>
    ompr::add_constraint(lam[p] <= 0, p = np + 1L) |>
    ompr::add_constraint(s[q] <= 0, q = nsq + 1L) |>
    ompr::add_constraint(a[i] <= 0, i = naq + 1L) |>
    ompr::add_constraint(w[m] <= 0, m = nwq + 1L) |>
    ompr::add_constraint(y[g] <= 0, g = ny + 1L) |>
    ompr::set_objective(
      ompr::sum_over(a_cost[i] * a[i], i = seq_len(naq)) +
        ompr::sum_over(w_cost[m] * w[m], m = seq_len(nwq)),
      sense = "min"
    )

  # z index of (tract c, team d) = (d - 1) * nu + c
  if (naq > 0L) {
    model <- model |>
      ompr::add_constraint(a[i] <= z[(pd[a_p[i]] - 1L) * nu + a_c[i]], i = seq_len(naq))
  }
  for (p in seq_len(np)) {
    model <- model |>
      ompr::add_constraint(ompr::sum_over(a[i], i = which(a_p == p)) == lam[p])
  }
  for (d in seq_len(nd)) {
    model <- model |>
      ompr::add_constraint(ompr::sum_over(z[c], c = (d - 1L) * nu + seq_len(nu)) == 1) |>
      ompr::add_constraint(ompr::sum_over(lam[p], p = ps_by_d[[d]]) <= k)
    if (length(ps_by_d[[d]]) + length(sq_by_d[[d]]) > 0L) {
      model <- model |>
        ompr::add_constraint(
          ompr::sum_over(pW[p] * lam[p], p = ps_by_d[[d]]) +
            ompr::sum_over(s_load[q] * s[q], q = sq_by_d[[d]]) >= lower - fixed_load[d]
        ) |>
        ompr::add_constraint(
          ompr::sum_over(pW[p] * lam[p], p = ps_by_d[[d]]) +
            ompr::sum_over(s_load[q] * s[q], q = sq_by_d[[d]]) <= upper - fixed_load[d]
        )
    }
  }
  for (u in match(free_ids, ids)) {
    model <- model |>
      ompr::add_constraint(
        ompr::sum_over(lam[mem_p[m]], m = which(mem_t == u)) +
          ompr::sum_over(s[q], q = which(sq_u == u)) == 1
      )
  }
  if (nwq > 0L) {
    nonfix <- which(!w_fixed)
    model <- model |>
      ompr::add_constraint(w[m] <= s[w_q[m]], m = seq_len(nwq))
    if (length(nonfix) > 0L) {
      model <- model |>
        ompr::add_constraint(w[nonfix[j]] <= y[w_y[j]], j = seq_along(nonfix))
    }
  }
  for (g in seq_len(ny)) {
    ps <- which(pd[mem_p] == y_key[g, 2] & mem_t == y_key[g, 1])
    model <- model |>
      ompr::add_constraint(y[g] == ompr::sum_over(lam[mem_p[m]], m = ps))
  }
  for (q in seq_len(nsq)) {
    model <- model |>
      ompr::add_constraint(ompr::sum_over(w[m], m = which(w_q == q)) == s[q])
  }

  solved <- .solve_model_with_status(
    model = model, solver = solver, max_time = max_time, rel_tol = rel_tol, verbose = verbose
  )
  if (!.is_solution_status(solved$status)) {
    msg <- tryCatch(
      solved$result$additional_solver_output$ROI$message$status_message,
      error = function(e) NULL
    )
    status <- if (solved$status == "error" && is.character(msg) &&
      any(grepl("infeasible", msg, ignore.case = TRUE))) "infeasible" else solved$status
    out <- finish_empty(status)
    out$diagnostics$optimality_gap <- solved$optimality_gap
    return(out)
  }
  res <- solved$result
  fill <- function(sol, n) {
    out <- numeric(n)
    out[sol[[2]]] <- sol$value
    out
  }
  lam_v <- fill(ompr::get_solution(res, lam[p]), np + 1L)[seq_len(np)]
  s_v <- fill(ompr::get_solution(res, s[q]), nsq + 1L)[seq_len(nsq)]
  z_v <- fill(ompr::get_solution(res, z[c]), nu * nd)
  w_v_sol <- fill(ompr::get_solution(res, w[m]), nwq + 1L)[seq_len(nwq)]
  .extract_assignment(
    solved = solved, start_time = start_time, parameters = parameters,
    cand = cand, ids = ids, team_ids = team_ids, load = load, dmat = dmat,
    pr = pr, pd = pd, sq_u = sq_u, sq_d = sq_d, w_q = w_q, w_v = w_v, w_cost = w_cost,
    lam_v = lam_v, s_v = s_v, z_v = z_v, w_sol = w_v_sol,
    fixed = fixed, fixed_load = fixed_load, lower = lower, upper = upper,
    warnings = cand_warnings
  )
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

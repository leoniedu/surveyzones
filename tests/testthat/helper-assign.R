toy_instance <- function(n = 6) {
  ids <- LETTERS[seq_len(n)]
  points <- sf::st_as_sf(
    tibble::tibble(
      point_id = tolower(paste0(ids, "1")),
      tract_id = ids,
      lon = -38.5 + 0.0028 * (seq_len(n) - 1),
      lat = -13
    ),
    coords = c("lon", "lat"), crs = 4326
  )
  idx <- stats::setNames(seq_len(n), ids)
  durations <- tidyr::expand_grid(origin_id = ids, destination_id = ids) |>
    dplyr::filter(origin_id != destination_id) |>
    dplyr::mutate(
      distance = unname(abs(idx[origin_id] - idx[destination_id])) +
        unname(ifelse(idx[origin_id] > idx[destination_id], 0.5, 0))
    )
  list(
    units = tibble::tibble(tract_id = ids, expected_service_time = as.numeric(seq_len(n))),
    points = points,
    teams = tibble::tibble(team_id = c("t1", "t2")),
    eligible = tidyr::expand_grid(tract_id = ids, team_id = c("t1", "t2")),
    durations = durations
  )
}

# Exhaustive optimum when every tract is mandatory (no loose/fixed tracts):
# every exact cover of the tracts by candidate roteiros, every team assignment.
brute_force_objective <- function(x, k, tol, radius_m = 500, max_tracts = 5L) {
  ids <- x$units$tract_id
  load <- stats::setNames(x$units$expected_service_time, ids)
  teams <- x$teams$team_id
  mean_load <- sum(load) / length(teams)
  cand <- surveyzones_roteiro_candidates(x$points, x$eligible, ids, radius_m, max_tracts)
  dur <- function(o, d) {
    if (o == d) return(0)
    x$durations$distance[x$durations$origin_id == o & x$durations$destination_id == d]
  }
  covers <- list()
  rec <- function(remaining, chosen) {
    if (length(remaining) == 0L) {
      covers[[length(covers) + 1L]] <<- chosen
      return(invisible())
    }
    first <- remaining[1]
    for (i in seq_len(nrow(cand))) {
      tr <- cand$tract_ids[[i]]
      if (first %in% tr && all(tr %in% remaining)) rec(setdiff(remaining, tr), c(chosen, i))
    }
  }
  rec(ids, integer())
  best <- Inf
  for (cv in covers) {
    opts <- lapply(cv, \(i) cand$teams[[i]])
    if (any(lengths(opts) == 0L)) next
    grid <- do.call(expand.grid, c(opts, stringsAsFactors = FALSE))
    for (g in seq_len(nrow(grid))) {
      tm <- unlist(grid[g, ], use.names = FALSE)
      wl <- vapply(teams, \(d) sum(vapply(cv[tm == d], \(i) sum(load[cand$tract_ids[[i]]]), numeric(1))), numeric(1))
      nr <- vapply(teams, \(d) sum(tm == d), numeric(1))
      if (any(nr > k) || any(wl < (1 - tol) * mean_load - 1e-9) || any(wl > (1 + tol) * mean_load + 1e-9)) next
      cost <- sum(vapply(teams, \(d) {
        mine <- cv[tm == d]
        if (length(mine) == 0L) return(0)
        min(vapply(ids, \(o) sum(vapply(mine, \(i) dur(o, cand$center_tract_id[i]), numeric(1))), numeric(1)))
      }, numeric(1)))
      best <- min(best, cost)
    }
  }
  best
}

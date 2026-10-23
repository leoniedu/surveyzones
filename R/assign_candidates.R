#' Candidate Walkable Roteiros for Team Assignment
#'
#' Builds candidate roteiros (sets of tracts walkable from a centre household)
#' for a team-assignment model. From each household taken as centre, a tract
#' is "near" when at least one of its households lies within `radius_m`
#' (straight line); the whole tract then joins. Every subset of near tracts
#' containing the centre's tract (up to `max_tracts` tracts) whose members
#' share at least one eligible team becomes a candidate. Every tract also gets
#' a singleton candidate. Candidates are deduplicated on the sorted tract set,
#' keeping the smallest `span_m`.
#'
#' @param points An sf object of POINT geometries (one row per household) with
#'   `tract_id` and `point_id` columns. Households without coordinates are
#'   simply absent.
#' @param eligible A data frame with columns `tract_id` and `team_id`.
#' @param tract_ids Character vector of all tract ids, including tracts that
#'   have no points.
#' @param radius_m Positive number. Maximum straight-line distance in metres
#'   from the centre household to the nearest household of a member tract.
#' @param max_tracts Positive integer. Maximum number of tracts per roteiro.
#'
#' @return A tibble ordered by `roteiro_id` with columns `roteiro_id`,
#'   `center_point_id` (`NA` for tracts without coordinates), `center_tract_id`,
#'   `tract_ids` (list of sorted character vectors), `teams` (list of teams
#'   eligible for all member tracts) and `span_m` (largest distance from the
#'   centre to the nearest household of each member tract; 0 for singletons).
#'   The `"warnings"` attribute is a character vector.
#' @export
surveyzones_roteiro_candidates <- function(points,
                                           eligible,
                                           tract_ids,
                                           radius_m = 500,
                                           max_tracts = 5L) {
  if (!inherits(points, "sf") || !all(c("tract_id", "point_id") %in% names(points))) {
    cli::cli_abort(
      "{.arg points} must be an sf object with {.field tract_id} and {.field point_id} columns."
    )
  }
  if (!all(as.character(sf::st_geometry_type(points)) == "POINT")) {
    cli::cli_abort("{.arg points} must contain only POINT geometries.")
  }
  if (!all(c("tract_id", "team_id") %in% names(eligible))) {
    cli::cli_abort(
      "{.arg eligible} must contain {.field tract_id} and {.field team_id} columns."
    )
  }
  if (!is.numeric(radius_m) || length(radius_m) != 1L || is.na(radius_m) || radius_m <= 0) {
    cli::cli_abort("{.arg radius_m} must be a single positive number.")
  }
  if (!is.numeric(max_tracts) || length(max_tracts) != 1L || is.na(max_tracts) ||
    max_tracts < 1 || max_tracts != round(max_tracts)) {
    cli::cli_abort("{.arg max_tracts} must be a whole number >= 1.")
  }
  tract_ids <- unique(as.character(tract_ids))
  pt_tract <- as.character(points$tract_id)
  pt_id <- as.character(points$point_id)
  if (anyNA(pt_tract)) {
    cli::cli_abort("{.field tract_id} in {.arg points} must not contain {.val NA}.")
  }
  if (anyNA(pt_id) || anyDuplicated(pt_id)) {
    cli::cli_abort("{.field point_id} in {.arg points} must be unique and not {.val NA}.")
  }
  unknown <- setdiff(unique(pt_tract), tract_ids)
  if (length(unknown) > 0L) {
    cli::cli_abort(
      "{.arg points} has tracts not in {.arg tract_ids}: {.val {unknown}}."
    )
  }
  max_tracts <- as.integer(max_tracts)

  teams_of <- split(as.character(eligible$team_id), as.character(eligible$tract_id))
  teams_of <- lapply(teams_of, \(x) sort(unique(x)))
  tract_teams <- \(id) teams_of[[id]] %||% character()

  best <- new.env(parent = emptyenv())
  record <- function(members, teams, span, center_pt, center_tract) {
    key <- paste(members, collapse = "+")
    cur <- best[[key]]
    if (is.null(cur) || span < cur$span ||
      (span == cur$span && !is.na(center_pt) &&
        (is.na(cur$center_point_id) || center_pt < cur$center_point_id))) {
      best[[key]] <- list(
        center_point_id = center_pt, center_tract_id = center_tract,
        tract_ids = members, teams = teams, span_m = span
      )
    }
  }

  warnings_out <- character()
  pts_by_tract <- split(seq_along(pt_tract), pt_tract)

  for (id in tract_ids) {
    idx <- pts_by_tract[[id]]
    if (is.null(idx)) {
      warnings_out <- c(warnings_out, paste0("tract ", id, " has no coordinates; singleton roteiro"))
      centre <- NA_character_
    } else {
      centre <- min(pt_id[idx])
    }
    if (length(tract_teams(id)) == 0L) {
      warnings_out <- c(warnings_out, paste0("tract ", id, " has no eligible team"))
    }
    record(id, tract_teams(id), 0, centre, id)
  }

  seen <- new.env(parent = emptyenv())
  geom <- sf::st_geometry(points)
  for (i in order(pt_id)) {
    own <- pt_tract[i]
    d <- as.numeric(sf::st_distance(geom[i], geom))
    tract_min <- vapply(pts_by_tract, \(idx) min(d[idx]), numeric(1))
    near <- sort(setdiff(names(tract_min)[tract_min <= radius_m], own))
    near <- near[lengths(lapply(near, \(t) intersect(tract_teams(own), tract_teams(t)))) > 0L]
    if (length(near) == 0L || max_tracts < 2L) next
    sig <- paste(own, paste(near, tract_min[near], sep = "=", collapse = "|"), sep = "#")
    if (!is.null(seen[[sig]])) next
    seen[[sig]] <- TRUE

    extend <- function(members, teams, start) {
      for (j in start:length(near)) {
        new_teams <- intersect(teams, tract_teams(near[j]))
        if (length(new_teams) == 0L) next
        m <- c(members, near[j])
        record(
          sort(c(own, m)), sort(new_teams),
          max(tract_min[m]), pt_id[i], own
        )
        if (length(m) + 1L < max_tracts && j < length(near)) {
          extend(m, new_teams, j + 1L)
        }
      }
    }
    extend(character(), tract_teams(own), 1L)
  }

  rows <- as.list(best)
  set_keys <- names(rows)
  rows <- unname(rows)
  out <- tibble::tibble(
    center_point_id = vapply(rows, \(r) r$center_point_id, character(1)),
    center_tract_id = vapply(rows, \(r) r$center_tract_id, character(1)),
    tract_ids = unname(lapply(rows, \(r) r$tract_ids)),
    teams = unname(lapply(rows, \(r) r$teams)),
    span_m = vapply(rows, \(r) r$span_m, numeric(1)),
    set_key = set_keys
  ) |>
    dplyr::arrange(.data$center_tract_id, .data$set_key) |>
    dplyr::mutate(roteiro_id = sprintf("R%04d", dplyr::row_number())) |>
    dplyr::select(
      "roteiro_id", "center_point_id", "center_tract_id",
      "tract_ids", "teams", "span_m"
    )

  attr(out, "warnings") <- warnings_out
  out
}

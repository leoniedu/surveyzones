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
  if (!is.numeric(max_tracts) || length(max_tracts) != 1L || is.na(max_tracts) || max_tracts < 1) {
    cli::cli_abort("{.arg max_tracts} must be a single positive integer.")
  }
  tract_ids <- unique(as.character(tract_ids))

  pt_tract <- as.character(points$tract_id)
  pt_id <- as.character(points$point_id)
  teams_of <- split(as.character(eligible$team_id), as.character(eligible$tract_id))
  teams_of <- lapply(teams_of, \(x) sort(unique(x)))
  get_teams <- \(ids) {
    sets <- lapply(ids, \(id) teams_of[[id]] %||% character())
    sort(Reduce(intersect, sets))
  }

  rows <- list()
  warnings_out <- character()

  if (nrow(points) > 0L) {
    dist_m <- surveyzones_engine_haversine(units = "m")(points, points)
    pts_by_tract <- split(seq_along(pt_tract), pt_tract)
    tract_names <- names(pts_by_tract)

    for (i in seq_along(pt_id)) {
      own <- pt_tract[i]
      tract_min <- vapply(
        pts_by_tract, \(idx) min(dist_m[i, idx]), numeric(1)
      )
      near <- setdiff(tract_names[tract_min <= radius_m], own)
      near <- near[vapply(
        near, \(t) length(get_teams(c(own, t))) > 0L, logical(1)
      )]
      for (k in 0:min(length(near), max_tracts - 1L)) {
        subsets <- if (k == 0L) list(character()) else utils::combn(near, k, simplify = FALSE)
        for (s in subsets) {
          members <- c(own, s)
          teams <- get_teams(members)
          if (length(teams) == 0L) next
          rows[[length(rows) + 1L]] <- tibble::tibble(
            center_point_id = pt_id[i],
            center_tract_id = own,
            tract_ids = list(sort(members)),
            teams = list(teams),
            span_m = if (k == 0L) 0 else max(tract_min[s])
          )
        }
      }
    }
  }

  has_points <- unique(pt_tract)
  for (id in setdiff(tract_ids, has_points)) {
    warnings_out <- c(warnings_out, paste0("tract ", id, " has no coordinates; singleton roteiro"))
    rows[[length(rows) + 1L]] <- tibble::tibble(
      center_point_id = NA_character_,
      center_tract_id = id,
      tract_ids = list(id),
      teams = list(get_teams(id)),
      span_m = 0
    )
  }

  out <- dplyr::bind_rows(rows)
  if (nrow(out) == 0L) {
    out <- tibble::tibble(
      center_point_id = character(), center_tract_id = character(),
      tract_ids = list(), teams = list(), span_m = numeric()
    )
  }
  out <- out |>
    dplyr::mutate(set_key = purrr::map_chr(.data$tract_ids, \(x) paste(x, collapse = "+"))) |>
    dplyr::arrange(.data$set_key, .data$span_m, .data$center_point_id) |>
    dplyr::distinct(.data$set_key, .keep_all = TRUE) |>
    dplyr::arrange(.data$center_tract_id, .data$set_key) |>
    dplyr::mutate(roteiro_id = sprintf("R%04d", dplyr::row_number())) |>
    dplyr::select(
      "roteiro_id", "center_point_id", "center_tract_id",
      "tract_ids", "teams", "span_m"
    )

  attr(out, "warnings") <- warnings_out
  out
}

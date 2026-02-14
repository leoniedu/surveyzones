#' Create a Haversine Distance Engine
#'
#' Returns a function that computes great-circle distances between
#' sf POINT objects using [sf::st_distance()].  This is the default
#' engine used by [surveyzones_compute_sparse_distances()] when no
#' engine is supplied.
#'
#' @param units Character scalar.
#'   One of `"km"` (default) or `"m"`.
#'
#' @return A function with signature
#'   `function(origins, destinations)` that returns an n x m numeric
#'   matrix of distances in the requested units.
#'
#' @export
surveyzones_engine_haversine <- function(units = c("km", "m")) {
  units <- match.arg(units)
  divisor <- if (units == "km") 1000 else 1

  function(origins, destinations) {
    dist_mat <- sf::st_distance(origins, destinations)
    # st_distance returns a units matrix; convert to plain numeric
    matrix(as.numeric(dist_mat) / divisor, nrow = nrow(origins))
  }
}


#' Create an OSRM Distance Engine
#'
#' Returns a function that computes road-network travel times or distances
#' using the OSRM routing service via [osrm::osrmTable()].
#'
#' @param measure Character scalar.  One of `"duration"` (default, in
#'   minutes) or `"distance"` (in meters).
#' @param osrm.server OSRM server URL.  Defaults to the `osrm` package
#'   option (the public demo server).
#' @param osrm.profile Routing profile.  Defaults to the `osrm` package
#'   option (typically `"car"`).
#'
#' @return A function with signature
#'   `function(origins, destinations)` that returns an n x m numeric
#'   matrix of travel times (minutes) or distances (meters).
#'
#' @details
#' The public OSRM demo server has rate limits and should not be used
#' for large-scale computation.
#' For production use, set up a local OSRM server and pass its URL
#' via the `osrm.server` argument.
#'
#' @export
surveyzones_engine_osrm <- function(measure = c("duration", "distance"),
                                    osrm.server = NULL,
                                    osrm.profile = NULL) {
  rlang::check_installed("osrm", reason = "to use the OSRM distance engine")
  measure <- match.arg(measure)

  function(origins, destinations) {
    opts <- options()
    on.exit(options(opts), add = TRUE)

    if (!is.null(osrm.server)) options(osrm.server = osrm.server)
    if (!is.null(osrm.profile)) options(osrm.profile = osrm.profile)

    result <- osrm::osrmTable(src = origins, dst = destinations,
                              measure = measure)

    if (measure == "duration") result$durations else result$distances
  }
}


#' Validate a Distance Engine
#'
#' Checks that a function satisfies the distance engine contract:
#' takes two sf POINT inputs and returns a numeric matrix with the
#' correct dimensions.
#'
#' @param engine A function to validate.
#' @param call Caller environment for error reporting.
#'
#' @return `engine`, invisibly (called for side effects).
#'
#' @keywords internal
validate_engine <- function(engine, call = rlang::caller_env()) {
  if (!is.function(engine)) {
    cli::cli_abort(
      "{.arg engine} must be a function, not {.obj_type_friendly {engine}}.",
      call = call
    )
  }

  # Check arity: must accept at least 2 arguments
  fmls <- formals(engine)
  if (length(fmls) < 2) {
    cli::cli_abort(
      c(
        "{.arg engine} must accept at least 2 arguments (origins, destinations).",
        "i" = "See {.fun surveyzones_engine_haversine} for the expected signature."
      ),
      call = call
    )
  }

  invisible(engine)
}

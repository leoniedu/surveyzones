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


#' Create a Google Maps Distance Matrix Engine
#'
#' Returns a function that computes travel times or distances using
#' the Google Maps Distance Matrix API via
#' [googleway::google_distance()].
#'
#' @param measure Character scalar.  One of `"duration"` (default, in
#'   minutes) or `"distance"` (in km).
#' @param mode Character scalar.  One of `"driving"` (default),
#'   `"walking"`, `"bicycling"`, or `"transit"`.
#' @param key Google Maps API key.  Defaults to the key set via
#'   [googleway::set_key()].
#' @param batch_size Integer scalar.  Maximum number of origins or
#'   destinations per API request.  The Google Distance Matrix API
#'   allows at most 25 origins or 25 destinations per call; results
#'   are stitched together transparently.
#'
#' @return A function with signature
#'   `function(origins, destinations)` that returns an n x m numeric
#'   matrix of travel times (minutes) or distances (km).
#'
#' @details
#' The Google Distance Matrix API is a paid service.  Each element
#' (origin-destination pair) is billed; see
#' <https://developers.google.com/maps/documentation/distance-matrix/usage-and-billing>
#' for current pricing.
#'
#' @export
surveyzones_engine_google <- function(measure = c("duration", "distance"),
                                      mode = c("driving", "walking",
                                               "bicycling", "transit"),
                                      key = NULL,
                                      batch_size = 25L) {
  rlang::check_installed("googleway",
                         reason = "to use the Google Maps distance engine")
  measure <- match.arg(measure)
  mode <- match.arg(mode)
  batch_size <- as.integer(batch_size)

  function(origins, destinations) {
    # Extract coordinates as lat/lon data.frames
    orig_coords <- sf::st_coordinates(origins)
    dest_coords <- sf::st_coordinates(destinations)
    orig_df <- data.frame(lat = orig_coords[, 2], lon = orig_coords[, 1])
    dest_df <- data.frame(lat = dest_coords[, 2], lon = dest_coords[, 1])

    n_orig <- nrow(orig_df)
    n_dest <- nrow(dest_df)

    # Batch origins and destinations to stay within API limits
    orig_batches <- split(seq_len(n_orig),
                          ceiling(seq_len(n_orig) / batch_size))
    dest_batches <- split(seq_len(n_dest),
                          ceiling(seq_len(n_dest) / batch_size))

    result_mat <- matrix(NA_real_, nrow = n_orig, ncol = n_dest)

    for (oi in orig_batches) {
      for (di in dest_batches) {
        api_args <- list(
          origins      = orig_df[oi, , drop = FALSE],
          destinations = dest_df[di, , drop = FALSE],
          mode         = mode
        )
        if (!is.null(key)) api_args$key <- key

        res <- do.call(googleway::google_distance, api_args)

        if (!identical(res$status, "OK")) {
          cli::cli_abort(
            "Google Distance Matrix API returned status {.val {res$status}}."
          )
        }

        # Parse the nested rows/elements structure
        for (i in seq_along(oi)) {
          elements <- res$rows$elements[[i]]
          if (measure == "duration") {
            values <- elements$duration$value  # seconds
            values <- values / 60              # -> minutes
          } else {
            values <- elements$distance$value  # meters
            values <- values / 1000            # -> km
          }
          # NA for non-OK element status
          statuses <- elements$status
          values[statuses != "OK"] <- NA_real_
          result_mat[oi[i], di] <- values
        }
      }
    }

    result_mat
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

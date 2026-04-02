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


#' Create a Google Routes API Distance Engine
#'
#' Returns a function that computes travel times or distances using
#' the Google Routes API `computeRouteMatrix` endpoint via
#' [httr2][httr2::httr2-package].
#'
#' @param measure Character scalar.  One of `"duration"` (default, in
#'   minutes) or `"distance"` (in km).
#' @param travel_mode Character scalar.  One of `"DRIVE"` (default),
#'   `"WALK"`, `"BICYCLE"`, or `"TRANSIT"`.
#' @param key Google Maps API key.  If `NULL`, reads from the
#'   environment variable `GOOGLE_MAPS_API_KEY`.
#' @param batch_size Integer scalar.  Maximum number of origins or
#'   destinations per API request.  The Routes API allows up to 25
#'   origins and 25 destinations per call (625 elements); results are
#'   stitched together transparently.
#'
#' @return A function with signature
#'   `function(origins, destinations)` that returns an n x m numeric
#'   matrix of travel times (minutes) or distances (km).
#'
#' @details
#' Uses the Routes API `computeRouteMatrix` endpoint
#' (<https://developers.google.com/maps/documentation/routes/compute_route_matrix>).
#' Each element (origin-destination pair) is billed; see
#' <https://developers.google.com/maps/documentation/routes/usage-and-billing>
#' for current pricing.
#'
#' Requires the **Routes API** to be enabled in your Google Cloud
#' project (not the legacy Distance Matrix API).
#'
#' @export
surveyzones_engine_google <- function(measure = c("duration", "distance"),
                                      travel_mode = c("DRIVE", "WALK",
                                                      "BICYCLE", "TRANSIT"),
                                      key = NULL,
                                      batch_size = 25L) {
  rlang::check_installed("httr2",
                         reason = "to use the Google Routes distance engine")
  rlang::check_installed("jsonlite",
                         reason = "to use the Google Routes distance engine")
  measure <- match.arg(measure)
  travel_mode <- match.arg(travel_mode)
  batch_size <- as.integer(batch_size)

  key <- key %||% Sys.getenv("GOOGLE_MAPS_API_KEY", unset = "")
  if (!nzchar(key)) {
    cli::cli_abort(c(
      "No Google Maps API key found.",
      "i" = "Set {.envvar GOOGLE_MAPS_API_KEY} or pass {.arg key} directly."
    ))
  }

  field_mask <- if (measure == "duration") {
    "originIndex,destinationIndex,duration"
  } else {
    "originIndex,destinationIndex,distanceMeters"
  }

  function(origins, destinations) {
    orig_coords <- sf::st_coordinates(origins)
    dest_coords <- sf::st_coordinates(destinations)

    n_orig <- nrow(orig_coords)
    n_dest <- nrow(dest_coords)

    orig_batches <- split(seq_len(n_orig),
                          ceiling(seq_len(n_orig) / batch_size))
    dest_batches <- split(seq_len(n_dest),
                          ceiling(seq_len(n_dest) / batch_size))

    result_mat <- matrix(NA_real_, nrow = n_orig, ncol = n_dest)

    for (oi in orig_batches) {
      for (di in dest_batches) {
        body <- list(
          origins = lapply(oi, function(i) {
            list(waypoint = list(location = list(latLng = list(
              latitude  = orig_coords[i, 2],
              longitude = orig_coords[i, 1]
            ))))
          }),
          destinations = lapply(di, function(j) {
            list(waypoint = list(location = list(latLng = list(
              latitude  = dest_coords[j, 2],
              longitude = dest_coords[j, 1]
            ))))
          }),
          travelMode = travel_mode
        )

        resp <- httr2::request(
          "https://routes.googleapis.com/distanceMatrix/v2:computeRouteMatrix"
        ) |>
          httr2::req_headers(
            `X-Goog-Api-Key`  = key,
            `X-Goog-FieldMask` = field_mask,
            `Content-Type`    = "application/json"
          ) |>
          httr2::req_body_json(body) |>
          httr2::req_error(is_error = function(r) FALSE) |>
          httr2::req_perform()

        if (httr2::resp_status(resp) != 200L) {
          cli::cli_abort(
            "Routes API returned HTTP {.val {httr2::resp_status(resp)}}: {httr2::resp_body_string(resp)}"
          )
        }

        # Response is a JSON array of route matrix elements
        elements <- jsonlite::fromJSON(httr2::resp_body_string(resp))

        for (row_idx in seq_len(nrow(elements))) {
          el <- elements[row_idx, ]
          # API returns 0-based indices
          oi_idx <- el$originIndex + 1L
          di_idx <- el$destinationIndex + 1L

          if (measure == "duration") {
            # duration is a string like "123s"
            secs <- as.numeric(sub("s$", "", el$duration))
            result_mat[oi[oi_idx], di[di_idx]] <- secs / 60
          } else {
            result_mat[oi[oi_idx], di[di_idx]] <- el$distanceMeters / 1000
          }
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

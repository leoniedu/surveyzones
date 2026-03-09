#' Compute All Pairwise Distances
#'
#' Computes distances between tract access points using a pluggable engine.
#' Distances are computed in chunks to avoid materialising a full N x N matrix.
#' Returns all non-self pairs; filtering by distance threshold should be done
#' by the caller (e.g., in [surveyzones_build_zones()]) to allow reuse with
#' different distance thresholds without recomputation.
#'
#' @param access_points An sf object with POINT geometries.
#'   Must contain a `tract_id` column.
#' @param engine A distance engine function (see
#'   [surveyzones_engine_haversine()]).  Default uses haversine in km.
#' @param chunk_size Integer.  Number of origins per chunk.
#'
#' @return A tibble with columns `origin_id`, `destination_id`,
#'   and `distance`.  All non-self pairs are included; self-pairs are
#'   excluded. For asymmetric engines (e.g., OSRM), both directional pairs
#'   are included (if present in the engine output).
#'
#' @export
surveyzones_compute_sparse_distances <- function(
    access_points,
    engine = surveyzones_engine_haversine(),
    chunk_size = 100L) {

  validate_engine(engine)
  validate_access_points(access_points)

  tract_ids <- access_points$tract_id
  n <- length(tract_ids)
  chunks <- split(seq_len(n), ceiling(seq_len(n) / chunk_size))

  results <- vector("list", length(chunks))

  cli::cli_progress_bar(
    "Computing distances",
    total = length(chunks),
    .envir = parent.frame()
  )

  for (i in seq_along(chunks)) {
    idx <- chunks[[i]]
    origins <- access_points[idx, ]
    # engine returns an n_chunk x N matrix
    dist_mat <- engine(origins, access_points)

    # Convert to long format (no filtering)
    chunk_dt <- .matrix_to_sparse_dt(
      dist_mat,
      origin_ids = tract_ids[idx],
      destination_ids = tract_ids,
      D_max = Inf
    )

    results[[i]] <- chunk_dt
    cli::cli_progress_update(.envir = parent.frame())
  }

  cli::cli_progress_done(.envir = parent.frame())

  dt <- dplyr::bind_rows(results)
  # Ensure no duplicates from overlapping chunks (shouldn't happen, but safe)
  dt |> dplyr::distinct(origin_id, destination_id, .keep_all = TRUE)
}


#' Wrap Pre-computed Distances
#'
#' Wraps an existing distance table into the sparse distance format
#' used by surveyzones, filtering to `D_max`.  Use this when you
#' already have distances from an external source (e.g. orce, a
#' database, or a pre-built dodgr matrix).
#'
#' @param distance_table A data.frame with columns `origin_id`,
#'   `destination_id`, and `distance`.
#' @param D_max Numeric scalar.  Maximum distance to retain.
#'
#' @return A tibble with the same schema, filtered and sorted.
#'
#' @export
surveyzones_precomputed_distances <- function(distance_table, D_max) {
  required <- c("origin_id", "destination_id", "distance")
  missing_cols <- setdiff(required, names(distance_table))
  if (length(missing_cols) > 0) {
    cli::cli_abort(
      "{.arg distance_table} is missing column{?s}: {.val {missing_cols}}."
    )
  }

  if (!is.numeric(D_max) || length(D_max) != 1 || D_max <= 0) {
    cli::cli_abort("{.arg D_max} must be a positive number.")
  }

  dt <- distance_table |>
    dplyr::as_tibble() |>
    dplyr::filter(distance <= D_max, origin_id != destination_id) |>
    dplyr::arrange(origin_id, destination_id)
  dt
}


#' Validate Access Points
#'
#' @param access_points Object to validate.
#' @param call Caller environment for error reporting.
#'
#' @return `access_points`, invisibly.
#' @keywords internal
validate_access_points <- function(access_points,
                                   call = rlang::caller_env()) {
  if (!inherits(access_points, "sf")) {
    cli::cli_abort(
      "{.arg access_points} must be an sf object.",
      call = call
    )
  }

  if (!"tract_id" %in% names(access_points)) {
    cli::cli_abort(
      "{.arg access_points} must contain a {.field tract_id} column.",
      call = call
    )
  }

  geom_types <- unique(as.character(sf::st_geometry_type(access_points)))
  if (!all(geom_types %in% "POINT")) {
    cli::cli_abort(
      "{.arg access_points} must contain only POINT geometries, not {.val {geom_types}}.",
      call = call
    )
  }

  if (anyDuplicated(access_points$tract_id)) {
    cli::cli_abort(
      "{.field tract_id} values in {.arg access_points} must be unique.",
      call = call
    )
  }

  invisible(access_points)
}


#' Convert Distance Matrix to Sparse Tibble
#'
#' @param mat Numeric matrix (n_origins x n_destinations).
#' @param origin_ids Character/integer vector of origin identifiers.
#' @param destination_ids Character/integer vector of destination identifiers.
#' @param D_max Maximum distance threshold.
#'
#' @return A tibble with columns `origin_id`, `destination_id`,
#'   `distance`.
#' @keywords internal
.matrix_to_sparse_dt <- function(mat, origin_ids, destination_ids, D_max) {
  # Find pairs within D_max (excluding self-pairs and NAs)
  within <- which(!is.na(mat) & mat <= D_max, arr.ind = TRUE)

  if (nrow(within) == 0L) {
    return(tibble::tibble(
      origin_id = character(0L),
      destination_id = character(0L),
      distance = numeric(0L)
    ))
  }

  oi <- as.character(origin_ids[within[, 1L]])
  di <- as.character(destination_ids[within[, 2L]])

  dt <- tibble::tibble(
    origin_id = oi,
    destination_id = di,
    distance = mat[within]
  )

  # Remove self-pairs
  dt |> dplyr::filter(origin_id != destination_id)
}


#' Complete a Sparse Distance Table with Haversine Fill-In
#'
#' For every directed pair of tract IDs in `access_points` not already
#' present in `sparse_distances`, computes haversine (great-circle) distance
#' in km and converts to the distance unit using a speed parameter:
#' `distance = haversine_km / speed_kmh * 60` (result in minutes).
#'
#' Existing pairs are preserved unchanged.  Use this before
#' [surveyzones_sequence()] when the
#' sparse distance table may be incomplete.
#'
#' @param sparse_distances A tibble with columns `origin_id`,
#'   `destination_id`, and `distance`.
#' @param access_points An sf object with POINT geometries and a
#'   `tract_id` column.
#' @param speed_kmh Numeric scalar.  Assumed travel speed in km/h for
#'   converting haversine distance to time.  Default `0.1` (intentionally
#'   harsh — missing pairs likely represent real barriers).
#'
#' @return A tibble with the same schema as `sparse_distances`, with
#'   missing pairs filled in.  Self-pairs are excluded.
#'
#' @export
surveyzones_complete_distances <- function(sparse_distances, access_points,
                                           speed_kmh = 0.1) {
  validate_access_points(access_points)

  if (!is.numeric(speed_kmh) || length(speed_kmh) != 1 || speed_kmh <= 0) {
    cli::cli_abort("{.arg speed_kmh} must be a positive number.")
  }

  tract_ids <- as.character(access_points$tract_id)
  n <- length(tract_ids)

  # All directed non-self pairs
  all_pairs <- tidyr::expand_grid(
    origin_id = tract_ids,
    destination_id = tract_ids
  ) |>
    dplyr::filter(origin_id != destination_id)

  # Find which pairs are missing
  existing <- sparse_distances |>
    dplyr::select(origin_id, destination_id) |>
    dplyr::mutate(
      origin_id = as.character(origin_id),
      destination_id = as.character(destination_id)
    )

  missing <- dplyr::anti_join(all_pairs, existing,
                               by = c("origin_id", "destination_id"))

  if (nrow(missing) == 0) {
    return(sparse_distances)
  }

  # Compute haversine for missing pairs
  idx_o <- match(missing$origin_id, tract_ids)
  idx_d <- match(missing$destination_id, tract_ids)

  haversine_km <- sf::st_distance(
    access_points[idx_o, ],
    access_points[idx_d, ],
    by_element = TRUE
  )
  haversine_km <- as.numeric(haversine_km) / 1000

  filled <- tibble::tibble(
    origin_id = missing$origin_id,
    destination_id = missing$destination_id,
    distance = haversine_km / speed_kmh * 60
  )

  cli::cli_alert_info(
    "Filled {nrow(filled)} missing distance pair{?s} using haversine (speed = {speed_kmh} km/h)."
  )

  dplyr::bind_rows(sparse_distances, filled)
}

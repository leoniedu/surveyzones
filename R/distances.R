#' Compute Sparse Pairwise Distances
#'
#' Computes distances between tract access points using a pluggable
#' engine, keeping only pairs within `D_max`.  Distances are computed
#' in chunks to avoid materialising a full N x N matrix.
#'
#' @param access_points An sf object with POINT geometries.
#'   Must contain a `tract_id` column.
#' @param D_max Numeric scalar.  Maximum distance to retain.
#' @param engine A distance engine function (see
#'   [surveyzones_engine_haversine()]).  Default uses haversine in km.
#' @param chunk_size Integer.  Number of origins per chunk.
#'
#' @return A `data.table` with columns `origin_id`, `destination_id`,
#'   and `travel_time`.  Only pairs with `travel_time <= D_max` are
#'   included; self-pairs are excluded.
#'
#' @export
surveyzones_compute_sparse_distances <- function(
    access_points,
    D_max,
    engine = surveyzones_engine_haversine(),
    chunk_size = 100L) {

  validate_engine(engine)
  validate_access_points(access_points)

  if (!is.numeric(D_max) || length(D_max) != 1 || D_max <= 0) {
    cli::cli_abort("{.arg D_max} must be a positive number.")
  }

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

    # Convert to long format and filter
    chunk_dt <- .matrix_to_sparse_dt(
      dist_mat,
      origin_ids = tract_ids[idx],
      destination_ids = tract_ids,
      D_max = D_max
    )

    results[[i]] <- chunk_dt
    cli::cli_progress_update(.envir = parent.frame())
  }

  cli::cli_progress_done(.envir = parent.frame())

  dt <- data.table::rbindlist(results)
  # Ensure no duplicates from overlapping chunks (shouldn't happen, but safe)
  data.table::setkey(dt, origin_id, destination_id)
  unique(dt)
}


#' Wrap Pre-computed Distances
#'
#' Wraps an existing distance table into the sparse distance format
#' used by surveyzones, filtering to `D_max`.  Use this when you
#' already have distances from an external source (e.g. orce, a
#' database, or a pre-built dodgr matrix).
#'
#' @param distance_table A data.frame with columns `origin_id`,
#'   `destination_id`, and `travel_time`.
#' @param D_max Numeric scalar.  Maximum distance to retain.
#'
#' @return A `data.table` with the same schema, filtered and keyed.
#'
#' @export
surveyzones_precomputed_distances <- function(distance_table, D_max) {
  required <- c("origin_id", "destination_id", "travel_time")
  missing_cols <- setdiff(required, names(distance_table))
  if (length(missing_cols) > 0) {
    cli::cli_abort(
      "{.arg distance_table} is missing column{?s}: {.val {missing_cols}}."
    )
  }

  if (!is.numeric(D_max) || length(D_max) != 1 || D_max <= 0) {
    cli::cli_abort("{.arg D_max} must be a positive number.")
  }

  dt <- data.table::as.data.table(distance_table)
  dt <- dt[travel_time <= D_max & origin_id != destination_id]
  data.table::setkey(dt, origin_id, destination_id)
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


#' Convert Distance Matrix to Sparse data.table
#'
#' @param mat Numeric matrix (n_origins x n_destinations).
#' @param origin_ids Character/integer vector of origin identifiers.
#' @param destination_ids Character/integer vector of destination identifiers.
#' @param D_max Maximum distance threshold.
#'
#' @return A `data.table` with columns `origin_id`, `destination_id`,
#'   `travel_time`.
#' @keywords internal
.matrix_to_sparse_dt <- function(mat, origin_ids, destination_ids, D_max) {
  # Find pairs within D_max (excluding self-pairs and NAs)
  within <- which(!is.na(mat) & mat <= D_max, arr.ind = TRUE)

  if (nrow(within) == 0L) {
    return(data.table::data.table(
      origin_id = character(0L),
      destination_id = character(0L),
      travel_time = numeric(0L)
    ))
  }

  oi <- as.character(origin_ids[within[, 1L]])
  di <- as.character(destination_ids[within[, 2L]])

  dt <- data.table::data.table(
    origin_id = oi,
    destination_id = di,
    travel_time = mat[within]
  )

  # Remove self-pairs
 dt[origin_id != destination_id]
}

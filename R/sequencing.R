#' Sequence All Zones in a Plan
#'
#' Computes a visit order for every zone using a TSP solver.
#'
#' @param plan A `surveyzones_plan` object.
#' @param sparse_distances Sparse distance table (as returned by
#'   [surveyzones_compute_sparse_distances()]).
#' @param method Character scalar.  TSP solving method passed to
#'   [TSP::solve_TSP()].  Default `"nn"` (nearest neighbour).
#'   Other useful options: `"repetitive_nn"`, `"nearest_insertion"`,
#'   `"cheapest_insertion"`, `"two_opt"`.
#'
#' @return The same `plan` with `plan$sequence` populated — a tibble
#'   with columns `zone_id`, `tract_id`, `visit_order`.
#'
#' @export
surveyzones_sequence <- function(plan, sparse_distances,
                                 method = "nn") {
  if (!inherits(plan, "surveyzones_plan")) {
    cli::cli_abort("{.arg plan} must be a {.cls surveyzones_plan} object.")
  }

  # Process each zone using purrr::map_df
  plan$sequence <- plan$assignments |>
    tidyr::nest(data = -zone_id) |>
    dplyr::mutate(
      ordered = purrr::map(
        data,
        \(z) surveyzones_sequence_zone(
          tract_ids = z$tract_id,
          sparse_distances = sparse_distances,
          method = method
        )
      )
    ) |>
    dplyr::select(-data) |>
    tidyr::unnest_longer(ordered, values_to = "tract_id") |>
    dplyr::group_by(zone_id) |>
    dplyr::mutate(visit_order = dplyr::row_number()) |>
    dplyr::ungroup()

  plan
}


#' Check if a Matrix is Symmetric
#'
#' @param mat Numeric matrix.
#' @param tol Tolerance for floating point comparison.
#'
#' @return Logical scalar.
#' @keywords internal
.is_symmetric_matrix <- function(mat, tol = 1e-8) {
  if (nrow(mat) != ncol(mat)) return(FALSE)
  all(is.na(mat) == is.na(t(mat))) &&
    isTRUE(all.equal(mat, t(mat), tolerance = tol, check.attributes = FALSE))
}


#' Build a Dense Distance Matrix for a Zone
#'
#' @param tract_ids Character vector of tract identifiers.
#' @param sparse_distances Sparse distance tibble.
#'
#' @return A numeric matrix with row/col names = tract_ids.
#'   Missing pairs filled with `Inf`. May be asymmetric.
#' @keywords internal
.build_zone_distance_matrix <- function(tract_ids, sparse_distances) {
  n <- length(tract_ids)
  mat <- matrix(Inf, nrow = n, ncol = n,
                dimnames = list(tract_ids, tract_ids))
  diag(mat) <- 0

  # Filter to within-zone pairs
  dt <- sparse_distances |>
    dplyr::filter(origin_id %in% tract_ids, destination_id %in% tract_ids)

  if (nrow(dt) > 0) {
    # Vectorized matrix assignment using cbind indexing
    mat[cbind(as.character(dt$origin_id), as.character(dt$destination_id))] <- dt$distance
  }

  mat
}


#' Sequence Tracts Within a Single Zone
#'
#' Orders tracts by solving a Travelling Salesman Problem using the
#' TSP package. Uses symmetric TSP for symmetric matrices, asymmetric
#' TSP (ATSP) for asymmetric distance matrices (e.g., from routing engines).
#'
#' @param tract_ids Character vector of tract identifiers in this zone.
#' @param sparse_distances Sparse distance table.
#' @param start Optional tract_id to start from.  Default: first element.
#' @param method Character scalar.  TSP solving method passed to
#'   [TSP::solve_TSP()].  Default `"nn"`.
#'
#' @return Character vector of tract_ids in visit order.
#'
#' @export
surveyzones_sequence_zone <- function(tract_ids, sparse_distances,
                                      start = NULL, method = "nn") {
  n <- length(tract_ids)
  if (n <= 2) return(tract_ids)

  # Build a small distance matrix for this zone
  dist_mat <- .build_zone_distance_matrix(tract_ids, sparse_distances)

  # Convert to TSP object (use ATSP for asymmetric matrices)
  if (.is_symmetric_matrix(dist_mat)) {
    tsp <- TSP::TSP(dist_mat, labels = tract_ids)
  } else {
    tsp <- TSP::ATSP(dist_mat, labels = tract_ids)
  }

  # Solve
  start_idx <- if (!is.null(start)) match(start, tract_ids) else 1L
  tour <- TSP::solve_TSP(tsp, method = method, start = start_idx)

  # Extract ordered tract_ids
  tract_ids[as.integer(tour)]
}

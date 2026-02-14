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

  zones <- split(plan$assignments, plan$assignments$zone_id)

  results <- vector("list", length(zones))

  for (i in seq_along(zones)) {
    z <- zones[[i]]
    zone_id <- z$zone_id[1]
    tract_ids <- z$tract_id

    ordered <- surveyzones_sequence_zone(
      tract_ids = tract_ids,
      sparse_distances = sparse_distances,
      method = method
    )

    results[[i]] <- tibble::tibble(
      zone_id = zone_id,
      tract_id = ordered,
      visit_order = seq_along(ordered)
    )
  }

  plan$sequence <- do.call(rbind, results)
  plan
}


#' Sequence Tracts Within a Single Zone
#'
#' Orders tracts by solving a Travelling Salesman Problem using the
#' TSP package.
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

  # Convert to TSP object
  tsp <- TSP::TSP(dist_mat, labels = tract_ids)

  # Solve
  start_idx <- if (!is.null(start)) match(start, tract_ids) else 1L
  tour <- TSP::solve_TSP(tsp, method = method, start = start_idx)

  # Extract ordered tract_ids
  tract_ids[as.integer(tour)]
}


#' Build a Dense Distance Matrix for a Zone
#'
#' @param tract_ids Character vector of tract identifiers.
#' @param sparse_distances Sparse distance data.table.
#'
#' @return A symmetric numeric matrix with row/col names = tract_ids.
#'   Missing pairs filled with `Inf`.
#' @keywords internal
.build_zone_distance_matrix <- function(tract_ids, sparse_distances) {
  n <- length(tract_ids)
  mat <- matrix(Inf, nrow = n, ncol = n,
                dimnames = list(tract_ids, tract_ids))
  diag(mat) <- 0

  # Filter to within-zone pairs
  dt <- sparse_distances[
    origin_id %in% tract_ids & destination_id %in% tract_ids
  ]

  if (nrow(dt) > 0) {
    for (r in seq_len(nrow(dt))) {
      oi <- as.character(dt$origin_id[r])
      di <- as.character(dt$destination_id[r])
      mat[oi, di] <- dt$travel_time[r]
      mat[di, oi] <- dt$travel_time[r]  # symmetrise
    }
  }

  mat
}

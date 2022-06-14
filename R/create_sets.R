#' Create Riemann resamples
#'
#' @param data sf object
#' @param min_size numeric (min polygon radius)
#' @param max_size numeric (max polygon radius)
#' @param step numeric (interval for cellsize steps)
#' @param ... passed to [sf::st_make_grid()]
#'
#' @export
make_riab_sets <- function(data, min_size, max_size, step, ...) {
  rlang::check_dots_used()

  dots <- rlang::list2(...)
  if (!is.null(dots$cellsize) || !is.null(dots$n)) {
    rlang::abort(
      "`cellsize` and `n` must not be arguments in `...`"
    )
  }

  if (!"sf" %in% class(data)) {
    rlang::abort(
      c(
        "`make_riab_sets()` currently only supports `sf` objects.",
        i = "Try converting `data` to an `sf` object via `sf::st_as_sf()`."
      )
    )
  }

  if (sf::st_crs(data) == sf::NA_crs_) {
    rlang::abort(
      c(
        "`make_riab_sets()` requires your data to have an appropriate coordinate reference system (CRS).",
        i = "Try setting a CRS using `sf::st_set_crs()`."
      )
    )
  }

  if (sf::st_is_longlat(data) && !sf::sf_use_s2()) {
    rlang::abort(
      c(
        "`make_riab_sets()` can only process geographic coordinates when using the s2 geometry library",
        "i" = "Reproject your data into a projected coordinate reference system using `sf::st_transform()`",
        "i" = "Or install the `s2` package and enable it using `sf::sf_use_s2(TRUE)`"
      )
    )
  }

  centroids <- sf::st_centroid(sf::st_geometry(data))

  grid_box <- sf::st_bbox(data)
  if (sf::st_is_longlat(data)) {
    # cf https://github.com/ropensci/stplanr/pull/467
    # basically: spherical geometry means sometimes the straight line of the
    # grid will exclude points within the bounding box
    #
    # so here we'll expand our boundary by a small bit in order to always contain our
    # points within the grid
    grid_box <- expand_grid(grid_box)
  }

  split_objs <- make_riemann_sets(data, centroids, grid_box, min_size, max_size, step, ...)

  cv_att <- list(
    min_size = min_size,
    max_size = max_size,
    step = step,
    ...
  )

  new_rset(
    splits = split_objs$splits,
    ids = split_objs[, grepl("^id", names(split_objs))],
    attrib = cv_att,
    subclass = c("riab_set", "spatial_rset", "rset")
  )

}

make_riemann_sets <- function(data, centroids, grid_box, min_size, max_size, step, ...) {
  n <- nrow(data)
  cellsizes <- seq(min_size, max_size, by = step)

  sgbp_sets <- NULL
  for (i in seq_along(cellsizes)) {
    # This is going to grow this tibble in a nasty way
    sgbp_sets <- update_sgbp_sets(centroids, sgbp_sets, grid_box, cellsizes[[i]], ...)
  }

  sgbp_sets$sgbp <- lapply(sgbp_sets$sgbp, default_complement, n = n)

  split_objs <- purrr::map(
    sgbp_sets$sgbp,
    make_splits,
    data = data,
    class = c("riab_splits", "spatial_rsplit")
  )

  tibble::tibble(
    splits = split_objs,
    id_fold = names0(length(split_objs), "Fold"),
    purrr::map_dfr(sgbp_sets[c("id_min_size", "id_max_size")], as.character)
  )

}

update_sgbp_sets <- function(centroids, sgbp_sets, grid_box, cellsize, ...) {
  grid <- sf::st_make_grid(grid_box, cellsize = cellsize, ...)

  new_sets <- setdiff(
    unique(row_ids_intersecting_fold_blocks(grid, centroids)),
    list(integer(0))
  )

  if (!is.null(sgbp_sets)) {
    sgbp_sets$id_max_size[sgbp_sets$sgbp %in% new_sets] <- cellsize
    new_sets <- new_sets[!(new_sets %in% sgbp_sets$sgbp)]
  }

  new_sets <- tibble::tibble(
    sgbp = new_sets,
    id_min_size = cellsize,
    id_max_size = cellsize
  )

  dplyr::bind_rows(
    sgbp_sets,
    new_sets
  )
}



expand_grid <- function(grid_box, expansion = 0.00001) {
  grid_box[1] <- grid_box[1] - abs(grid_box[1] * expansion)
  grid_box[2] <- grid_box[2] - abs(grid_box[2] * expansion)
  grid_box[3] <- grid_box[3] + abs(grid_box[3] * expansion)
  grid_box[4] <- grid_box[4] + abs(grid_box[4] * expansion)
  grid_box
}


row_ids_intersecting_fold_blocks <- function(grid_blocks, data) {
  # grid_blocks is a list of sgbp lists (?sf::sgbp)
  #
  # The first map() here iterates through the meta-list,
  # and the second checks each element of the relevant sgbp list
  # to see if it is integer(0) (no intersections) or not
  #
  # Each sgbp sub-list is nrow(data) elements long, so this which()
  # returns the list indices which are not empty, which is equivalent
  # to the row numbers that intersect with blocks in the fold
  purrr::map(
    grid_blocks,
    function(blocks) {
      which(
        purrr::map_lgl(
          sf::st_intersects(data, blocks),
          sgbp_is_not_empty
        )
      )
    }
  )
}

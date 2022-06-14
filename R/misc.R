# Check sparse geometry binary predicate for empty elements
# See ?sf::sgbp for more information on the data structure
sgbp_is_not_empty <- function(x) !identical(x, integer(0))


`%?%` <- function(lhs, rhs) {
  if (!rlang::is_list(rhs, 2)) {
    rlang::abort("`rhs` must be a list of length 2", call = rlang::caller_env())
  }

  if (is.null(lhs) ||
      is.na(lhs) ||
      is.nan(lhs) ||
      is.infinite(lhs) ||
      (is.logical(lhs) && !lhs)) {
    rhs[[1]]
  } else {
    rhs[[2]]
  }
}


## Keep synced with rsample

names0 <- function(num, prefix = "x") {
  if (num == 0L) {
    return(character())
  }
  ind <- format(1:num)
  ind <- gsub(" ", "0", ind)
  paste0(prefix, ind)
}

default_complement <- function(ind, n) {
  list(
    analysis = setdiff(1:n, ind),
    assessment = unique(ind)
  )
}

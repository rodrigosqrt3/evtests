#' p-value from a quantile table of a limiting distribution
#'
#' `q` is a vector of quantiles at `limit_tables$probs`. Values beyond the
#' tabulated range are reported at its edges: p-values below
#' `1 - max(probs)` are returned as `1 - max(probs)`.
#' @noRd
table_pvalue <- function(stat, q) {
  probs <- limit_tables$probs
  f <- stats::approx(q, probs, xout = stat, yleft = 0, yright = max(probs),
                     ties = "ordered")$y
  max(1 - f, 1 - max(probs))
}

#' Quantile vector at a value of gamma, by linear interpolation between the
#' columns of a table (values outside the grid take the nearest column)
#' @noRd
table_at <- function(tab, g) {
  grid <- tab$gamma
  o <- order(grid)
  grid <- grid[o]
  q <- tab$q[, o, drop = FALSE]
  g <- min(max(g, grid[1L]), grid[length(grid)])
  j <- findInterval(g, grid, rightmost.closed = TRUE)  # in 1, ..., length - 1
  w <-(g - grid[j]) / (grid[j + 1L] - grid[j])
  (1 - w) * q[, j] + w * q[, j + 1L]
}

#' Construct a test path over k
#'
#' @param df Data frame with columns `k`, `statistic` and `p_value`, one row
#'   per value of k (other columns allowed).
#' @param method Name of the test.
#' @param null Description of the null hypothesis.
#' @param alternative The alternative.
#' @noRd
new_evt_test <- function(df, method, null, alternative) {
  rownames(df) <- NULL
  structure(df, class = c("evt_test", "data.frame"), method = method,
            null = null, alternative = alternative)
}

#' Validate the alternative of a test of the sign of gamma
#' @noRd
evt_alternative <- function(alternative) {
  match.arg(alternative, c("two.sided", "positive", "negative"))
}

#' Combine one-sided p-values according to the alternative
#'
#' @param p_pos p-value against `gamma > 0`.
#' @param p_neg p-value against `gamma < 0`.
#' @noRd
evt_pick <- function(p_pos, p_neg, alternative) {
  switch(alternative,
         positive = p_pos,
         negative = p_neg,
         two.sided = pmin(1, 2 * pmin(p_pos, p_neg)))
}

#' Print a test path
#'
#' @param x An `evt_test` object.
#' @param digits Number of significant digits.
#' @param ... Passed to the data frame print method.
#' @return `x`, invisibly.
#' @export
print.evt_test <- function(x, digits = 4L, ...) {
  cat("<", attr(x, "method"), ">\n", sep = "")
  cat("  null: ", attr(x, "null"), "; alternative: ",
      switch(attr(x, "alternative"), two.sided = "gamma != 0",
             positive = "gamma > 0", negative = "gamma < 0"), "\n", sep = "")
  print(as.data.frame(x), digits = digits, row.names = FALSE, ...)
  invisible(x)
}

#' Plot a test path against k
#'
#' Plots the p-value (default) or the statistic against `k`, with a
#' horizontal line at `alpha` for p-values.
#'
#' @param x An `evt_test` object with more than one value of k.
#' @param y Name of the column to plot, `"p_value"` by default.
#' @param alpha Level drawn as a horizontal line when plotting p-values.
#' @param ... Passed to [graphics::plot()].
#' @return `x`, invisibly.
#' @export
plot.evt_test <- function(x, y = "p_value", alpha = 0.05, ...) {
  if (!y %in% names(x)) stop("no column `", y, "` to plot.", call. = FALSE)
  graphics::plot(x$k, x[[y]], type = "l", xlab = "k", ylab = y,
                 main = attr(x, "method"), ...)
  if (grepl("^p_", y)) graphics::abline(h = alpha, lty = 2)
  invisible(x)
}

null_cache <- new.env(parent = emptyenv())

#' Simulated null distribution with a fixed seed
#'
#' Runs `gen(nsim)` (which must return a numeric vector of `nsim` simulated
#' statistics) with a fixed seed, caches the sorted result under `key`, and
#' restores the global random number stream.
#' @noRd
sim_null <- function(key, nsim, gen) {
  key <- paste(key, nsim)
  if (!is.null(null_cache[[key]])) return(null_cache[[key]])
  had_seed <- exists(".Random.seed", envir = globalenv(), inherits = FALSE)
  if (had_seed) old <- get(".Random.seed", envir = globalenv(), inherits = FALSE)
  on.exit({
    if (had_seed) {
      assign(".Random.seed", old, envir = globalenv())
    } else if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) {
      rm(".Random.seed", envir = globalenv())
    }
  }, add = TRUE)
  set.seed(20261120L)
  out <- sort(gen(nsim))
  null_cache[[key]] <- out
  out
}

#' Greenwood statistics sum(E^2) / sum(E)^2 of m exponentials, nsim times
#'
#' Generated in blocks of at most about 5e6 random numbers.
#' @noRd
greenwood_sim <- function(m, nsim) {
  block <- max(1L, floor(5e6 / m))
  out <- numeric(0)
  while (length(out) < nsim) {
    b <- min(block, nsim - length(out))
    e <- matrix(stats::rexp(b * m), b, m)
    out <- c(out, rowSums(e^2) / rowSums(e)^2)
  }
  out
}

#' Validate the number of simulations
#' @noRd
check_nsim <- function(nsim) {
  if (!is.numeric(nsim) || length(nsim) != 1L || is.na(nsim) ||
      nsim < 100 || nsim != round(nsim)) {
    stop("`nsim` must be a whole number of at least 100.", call. = FALSE)
  }
  invisible(nsim)
}

#' Excesses of the k largest observations over the (k + 1)-th largest
#' @noRd
evt_excess <- function(z, k) z[seq_len(k)] - z[k + 1L]

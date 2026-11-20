#' Hall & Wang (1999) minimum-distance end-point estimator
#'
#' Estimates a finite end-point from the extreme observations, which are
#' taken to follow exactly the limiting joint law of the extremes of Weibull
#' type, by the minimum-distance method of Hall & Wang (1999). The paper is
#' written for a lower end-point `theta` and the `m` smallest observations
#' `X[1] < ... < X[m]`: with
#' `Y[i] = i (log(X[i + 1] - theta) - log(X[i] - theta))`, `i = 1, ..., m - 1`
#' (equation 2.1), the estimate minimises Greenwood's statistic
#' `G(theta) = sum(Y^2) / sum(Y)^2` over `theta < X[1]` (Section 2.2). At the
#' true end-point the `Y[i]` are independent exponentials with a common mean,
#' so `G(theta_0)` is distributed as Greenwood's statistic of `m - 1`
#' exponentials, free of any parameter, and `{theta : G(theta) <= g}`, with
#' `g` the upper `alpha` point of that distribution, is a confidence region
#' with exact coverage `1 - alpha` (p. 185). An upper end-point is handled by
#' the transformation `x -> -x` (p. 179).
#'
#' **Convention for `k`.** As in the rest of the package, `k` is the number of
#' excesses: the estimator uses the `k + 1` most extreme observations (the
#' paper's `m = k + 1`).
#'
#' The shape `alpha` is estimated, as in the paper, from the generalized least
#' squares regression of `log(X[i] - theta_hat)` on the expected log-points of
#' the limiting process: `log Gamma[i]`, with mean `digamma(i)` and covariance
#' `trigamma(max(i, j))` (Weissman 1978), whose slope is `1 / alpha`.
#'
#' The null quantile is simulated (`nsim` draws, fixed seed). The estimator
#' does not exist when `G` has no minimum before the limit at `-Inf` (for a
#' lower end-point); the estimate is then `NA`. An interval limit that does
#' not exist is `-Inf` or `Inf`.
#'
#' @param x Numeric sample.
#' @param k Number of excesses, whole numbers of at least 3 (at least 4
#'   observations), or `NULL` for `min(n - 1, 11)`.
#' @param alpha Significance level, in (0, 1).
#' @param side `"upper"` (default) for an upper end-point, `"lower"` for a
#'   lower one, as in the paper.
#' @param nsim Number of simulations of the null distribution of `G`.
#'
#' @return An `evt_estimate` data frame with columns `k`, `endpoint`,
#'   `lower`, `upper` (the confidence interval), `shape` (the estimate of
#'   `alpha`), `statistic` (the minimum of `G`) and `critical`.
#'
#' @references
#' Hall, P., & Wang, J. Z. (1999). Estimating the end-point of a probability
#'   distribution using minimum-distance methods. Bernoulli, 5(1), 177-189.
#'   doi:10.2307/3318618
#'
#' @seealso [olympic_100m], the example of the paper.
#'
#' @examples
#' # Section 3 of the paper: 12 best 100 m times, lower end-point, 90%
#' hall_wang1999(olympic_100m$time, k = 11, alpha = 0.1, side = "lower",
#'               nsim = 2e4)
#'
#' @export
hall_wang1999 <- function(x, k = NULL, alpha = 0.05,
                          side = c("upper", "lower"), nsim = 1e5) {
  side <- match.arg(side)
  check_alpha(alpha)
  check_nsim(nsim)
  label <- "Hall & Wang (1999)"
  sgn <- if (side == "upper") 1 else -1
  z <- evt_sorted(sgn * x, label, min_n = 4L)
  n <- length(z)
  if (is.null(k)) k <- min(n - 1L, 11L)
  k <- evt_k(k, 3L, n - 1L, label)
  rows <- lapply(k, function(kk) {
    top <- z[seq_len(kk + 1L)]
    if (!(top[1L] > top[kk + 1L])) {
      stop("the k + 1 most extreme observations must not all be equal.",
           call. = FALSE)
    }
    g <- top[1L] - top                       # distances below the maximum
    nul <- sim_null(paste("greenwood", kk), nsim,
                    function(ns) greenwood_sim(kk, ns))
    crit <- stats::quantile(nul, 1 - alpha, names = FALSE)
    fit <- hw_fit(g)
    ci <- hw_interval(g, crit)
    shape <- if (is.na(fit$s)) NA_real_ else hw_shape(fit$s + g)
    est <- top[1L] + fit$s
    lim <- top[1L] + ci
    if (sgn < 0) {
      est <- -est
      lim <- rev(-lim)
    }
    data.frame(k = kk, endpoint = est, lower = lim[1L], upper = lim[2L],
               shape = shape, statistic = fit$stat, critical = crit)
  })
  new_evt_estimate(do.call(rbind, rows), label, param = "endpoint")
}

#' Greenwood's statistic G as a function of the distance s of the candidate
#' end-point beyond the most extreme observation (g: distances below it)
#' @noRd
hw_g <- function(s, g) {
  k <- length(g)
  logs <- log(outer(s, g, "+"))
  y <- logs[, -1L, drop = FALSE] - logs[, -k, drop = FALSE]
  y <- sweep(y, 2L, seq_len(k - 1L), "*")
  rowSums(y^2) / rowSums(y)^2
}

#' Limit of G as the end-point goes to infinity
#' @noRd
hw_g_inf <- function(g) {
  y <- seq_len(length(g) - 1L) * diff(g)
  sum(y^2) / sum(y)^2
}

#' Grid of distances, relative to the span of the observations
#' @noRd
hw_grid <- function(g) g[length(g)] * exp(seq(log(1e-6), log(1e6),
                                              length.out = 600L))

#' Minimum of G
#' @noRd
hw_fit <- function(g) {
  grid <- hw_grid(g)
  stat <- hw_g(grid, g)
  inf <- hw_g_inf(g)
  i <- which.min(stat)
  if (i > 1L && i < length(grid) && stat[i] < inf) {
    o <- stats::optimize(function(ls) hw_g(exp(ls), g),
                         log(c(grid[i - 1L], grid[i + 1L])), tol = 1e-12)
    list(s = exp(o$minimum), stat = o$objective)
  } else {
    list(s = NA_real_, stat = min(stat[i], inf))
  }
}

#' Confidence region {s : G(s) <= crit} as distances beyond the maximum
#' @noRd
hw_interval <- function(g, crit) {
  grid <- hw_grid(g)
  stat <- hw_g(grid, g)
  ok <- which(stat <= crit)
  inf_ok <- hw_g_inf(g) <= crit
  if (!length(ok) && !inf_ok) return(c(NA_real_, NA_real_))
  if (!length(ok)) return(c(grid[length(grid)], Inf))
  f <- function(ls) hw_g(exp(ls), g) - crit
  first <- ok[1L]
  last <- ok[length(ok)]
  lower <- if (first == 1L) 0 else
    exp(stats::uniroot(f, log(grid[c(first - 1L, first)]), tol = 1e-12)$root)
  upper <- if (last == length(grid)) Inf else
    exp(stats::uniroot(f, log(grid[c(last, last + 1L)]), tol = 1e-12)$root)
  c(lower, upper)
}

#' alpha from the GLS slope of log distances on E(log Gamma_i)
#' @noRd
hw_shape <- function(d) {
  m <- length(d)
  zz <- log(rev(sort(d, decreasing = TRUE)))   # smallest distance first
  ev <- digamma(seq_len(m))
  v <- trigamma(outer(seq_len(m), seq_len(m), pmax))
  q <- solve(v)
  xm <- cbind(ev, 1)
  beta <- solve(t(xm) %*% q %*% xm, t(xm) %*% q %*% zz)
  1 / beta[1L]
}

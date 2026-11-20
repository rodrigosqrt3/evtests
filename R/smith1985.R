#' Maximum likelihood end-point with a profile-likelihood interval (Smith
#' 1985)
#'
#' Fits the generalized Pareto distribution with a negative index to the `k`
#' excesses over `u = X[n-k,n]`, written in terms of the end-point `theta`:
#' `y / (theta - u)` follows a Beta(1, `a`) distribution, with shape
#' `a = -1 / gamma`. Maximising over `a` for fixed `theta` gives the profile
#' log-likelihood
#' `k log(a) - k log(theta - u) - k - sum(log((theta - X[i]) / (theta - u)))`
#' with `a = -k / sum(log((theta - X[i]) / (theta - u)))`, which is maximised
#' over `theta`. The interval is `{theta : 2 (l(theta_hat) - l(theta)) <=
#' qchisq(1 - alpha, 1)}`. This is the likelihood treatment of non-regular
#' end-point problems of Smith (1985); the estimate is the same as the
#' end-point `u + alpha / k` of the maximum likelihood fit of
#' [grimshaw1993()] when that fit has `k > 0` (Grimshaw's `k`).
#'
#' **Regularity.** The likelihood is unbounded as `theta` approaches the
#' sample maximum, so the estimate is an interior local maximum. It exists
#' and is asymptotically normal only when `a > 2` (`gamma > -1/2`; Smith
#' 1985, Theorems 1 and 3); for `1 < a <= 2` the problem is non-regular, and
#' for `a <= 1` there is no consistent maximum likelihood estimator. When the
#' fitted shape is not larger than 2 the function warns that the chi-squared
#' calibration of the interval is not justified; the interval of
#' [hall_wang1999()] is exact for any shape. The estimate is `NA` when the
#' profile likelihood has no interior maximum above its limit at infinity.
#' The paper has no numerical example.
#'
#' @param x Numeric sample.
#' @param k Number of excesses, a whole number of at least 4 and at most
#'   `n - 1`, or `NULL` for `min(n - 1, 10)`.
#' @param alpha Significance level, in (0, 1).
#'
#' @return An `evt_estimate` data frame with one row and columns `k`,
#'   `endpoint`, `lower`, `upper`, `shape` (`a`), `gamma` (`-1 / a`),
#'   `loglik` and `regular` (`TRUE` when the shape exceeds 2).
#'
#' @references
#' Smith, R. L. (1985). Maximum likelihood estimation in a class of
#'   nonregular cases. Biometrika, 72(1), 67-90. doi:10.1093/biomet/72.1.67
#'
#' @examples
#' set.seed(1)
#' x <- 1 - stats::runif(500)^(1 / 3)     # end-point 1, shape 3
#' smith1985(x, k = 50)
#'
#' @export
smith1985 <- function(x, k = NULL, alpha = 0.05) {
  check_alpha(alpha)
  label <- "Smith (1985)"
  z <- evt_sorted(x, label, min_n = 5L)
  n <- length(z)
  if (is.null(k)) k <- min(n - 1L, 10L)
  if (length(k) != 1L) stop("`k` must be a single number.", call. = FALSE)
  k <- evt_k(k, 4L, n - 1L, label)
  u <- z[k + 1L]
  top <- z[seq_len(k)]
  if (!(top[1L] > u)) {
    stop("the k + 1 largest observations must not all be equal.",
         call. = FALSE)
  }
  g <- top[1L] - u
  fit <- smith_fit(top, u, g)
  crit <- stats::qchisq(1 - alpha, 1) / 2
  ci <- if (is.na(fit$s)) c(NA_real_, NA_real_) else
    smith_interval(top, u, g, fit, crit)
  if (!is.na(fit$shape) && fit$shape <= 2) {
    warning("the fitted shape (", format(fit$shape, digits = 3), ") is not ",
            "larger than 2: the problem is non-regular (Smith 1985), and the ",
            "chi-squared calibration of the interval is not justified.",
            call. = FALSE)
  }
  out <- data.frame(k = k,
                    endpoint = top[1L] + fit$s,
                    lower = top[1L] + ci[1L], upper = top[1L] + ci[2L],
                    shape = fit$shape, gamma = -1 / fit$shape,
                    loglik = fit$ll, regular = isTRUE(fit$shape > 2))
  new_evt_estimate(out, label, param = "endpoint")
}

#' Profile log-likelihood of the end-point at distance s beyond the maximum
#' @noRd
smith_profile <- function(s, x, u) {
  k <- length(x)
  r <- (x[1L] - x + s) / (x[1L] - u + s)
  slr <- sum(log(r))
  a <- -k / slr
  k * log(a) - k * log(x[1L] - u + s) - k - slr
}

#' Limit of the profile at infinity (exponential excesses)
#' @noRd
smith_profile_inf <- function(x, u) -length(x) * log(mean(x - u)) - length(x)

#' Interior local maximum of the profile likelihood
#' @noRd
smith_fit <- function(x, u, g) {
  k <- length(x)
  grid <- g * exp(seq(log(1e-4), log(1e5), length.out = 600L))
  ll <- vapply(grid, smith_profile, numeric(1), x = x, u = u)
  m <- length(ll)
  i <- which(ll[-c(1L, m)] > ll[-c(m - 1L, m)] & ll[-c(1L, m)] > ll[-c(1L, 2L)]) +
    1L
  none <- list(s = NA_real_, shape = NA_real_, ll = NA_real_)
  if (!length(i)) return(none)
  i <- i[which.max(ll[i])]
  o <- stats::optimize(function(ls) smith_profile(exp(ls), x, u),
                       log(grid[c(i - 1L, i + 1L)]), maximum = TRUE,
                       tol = 1e-12)
  s <- exp(o$maximum)
  if (!(o$objective > smith_profile_inf(x, u))) return(none)
  r <- (x[1L] - x + s) / (x[1L] - u + s)
  list(s = s, shape = -k / sum(log(r)), ll = o$objective)
}

#' Profile-likelihood interval around the interior maximum
#' @noRd
smith_interval <- function(x, u, g, fit, crit) {
  target <- fit$ll - crit
  f <- function(ls) smith_profile(exp(ls), x, u) - target
  up <- fit$s * exp(seq(0, log(1e5 * g / fit$s), length.out = 400L))
  fu <- vapply(log(up), f, numeric(1))
  below <- which(fu < 0)
  upper <- if (length(below)) {
    j <- below[1L]
    exp(stats::uniroot(f, log(up[c(j - 1L, j)]), tol = 1e-12)$root)
  } else Inf
  lo <- fit$s * exp(seq(0, log(1e-6), length.out = 400L))
  fl <- vapply(log(lo), f, numeric(1))
  belowl <- which(fl < 0)
  lower <- if (length(belowl)) {
    j <- belowl[1L]
    exp(stats::uniroot(f, log(lo[c(j, j - 1L)]), tol = 1e-12)$root)
  } else 0
  c(lower, upper)
}

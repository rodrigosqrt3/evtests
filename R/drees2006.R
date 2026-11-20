#' Test of the extreme value condition of Drees, de Haan & Li (2006)
#'
#' Tests whether the distribution is in the domain of attraction of an
#' extreme value distribution (`H0`) by comparing the tail empirical
#' distribution function of the `k` largest observations with the fitted
#' generalized Pareto approximation:
#' `k T = k integral_0^1 ((n/k) Fbar_n(a (x^(-gamma) - 1) / gamma + b) -
#' x)^2 x^(eta - 2) dx`, with `b = X[n-k,n]` and `(gamma, a)` the maximum
#' likelihood estimates of the generalized Pareto distribution fitted to the
#' `k` excesses over `b` (Example 2.1; computed with [grimshaw1993()]), and
#' `eta = 1`. It works on the original scale (no logarithms) and is location
#' and scale invariant. The integral is evaluated exactly: with
#' `v[i] = (1 + gamma y[i] / a)^(-1 / gamma)` for the excesses `y[i]`, the
#' integrand is `(G(x) - x)^2 / x` with `G` the empirical distribution
#' function of the `v[i]`.
#'
#' **Null distribution.** Under `H0`, `gamma > -1/2` and `sqrt(k) A(n/k) -> 0`
#' (Theorem 2.2), `k T` converges to a functional of a Brownian motion whose
#' distribution depends on `gamma`. The paper rejects when `k T` exceeds the
#' quantile at the estimated `gamma`, interpolated linearly in its Table 1,
#' with the column `gamma = 2` for larger estimates, and does not apply the
#' test when the estimate is below -0.499 or does not exist. The quantiles
#' used here come from a simulation of the limiting variable (50000
#' replicates per value of `gamma` on a grid from -0.49 to 2, Brownian motion
#' on 4000 points), which reproduces Table 1 of the paper (see the
#' reproduction vignette for one misprinted cell); p-values below 5e-4 are
#' reported as 5e-4. The paper notes that the size grows when `k` exceeds the
#' value that is optimal for estimation, so that the test is a tool to find
#' the largest `k` for which the generalized Pareto approximation holds.
#'
#' @param x Numeric sample.
#' @param k Number of excesses: whole numbers from 2 to `n - 1`, or `NULL`
#'   for all of them.
#'
#' @return An `evt_test` data frame with columns `k`, `statistic` (`k T`),
#'   `p_value` (`NA` when the maximum likelihood estimate does not exist or
#'   is below -0.499) and `gamma`.
#'
#' @references
#' Drees, H., de Haan, L., & Li, D. (2006). Approximations to the tail
#'   empirical distribution function with application to testing extreme
#'   value conditions. Journal of Statistical Planning and Inference, 136(10),
#'   3498-3538. doi:10.1016/j.jspi.2005.02.017
#'
#' @examples
#' set.seed(1)
#' drees2006(stats::rnorm(1000), k = c(50, 100, 200))
#'
#' @export
drees2006 <- function(x, k = NULL) {
  label <- "Drees, de Haan & Li (2006)"
  z <- evt_sorted(x, label, min_n = 3L)
  n <- length(z)
  k <- evt_k(k, 2L, n - 1L, label)
  rows <- lapply(k, function(kk) {
    y <- evt_excess(z, kk)
    fit <- if (any(y > 0)) grimshaw1993(y) else list(exists = FALSE)
    if (!fit$exists) {
      return(data.frame(k = kk, statistic = NA_real_, p_value = NA_real_,
                        gamma = NA_real_))
    }
    s <- drees2006_stat(y, fit$gamma, fit$scale)
    p <- if (fit$gamma < -0.499) NA_real_ else
      table_pvalue(s, table_at(limit_tables$T, fit$gamma))
    data.frame(k = kk, statistic = s, p_value = p, gamma = fit$gamma)
  })
  new_evt_test(do.call(rbind, rows), paste0(label, ", eta = 1"),
               "extreme value condition", "two.sided")
}

#' k T of Drees, de Haan & Li (2006) with eta = 1, exact integral
#' @noRd
drees2006_stat <- function(y, g, a) {
  k <- length(y)
  v <- if (abs(g) < 1e-12) exp(-y / a) else (1 + g * y / a)^(-1 / g)
  v <- sort(pmin(pmax(v, 0), 1))
  # integral of (c - x)^2 / x over pieces where G = c = j / k; on (0, v[1])
  # G = 0 and the integrand is x
  lo <- c(0, v)
  hi <- c(v, 1)
  cc <- (0:k) / k
  prim <- function(x, c) c^2 * log(x) - 2 * c * x + x^2 / 2
  tot <- hi[1L]^2 / 2
  for (j in 2:(k + 1L)) {
    if (hi[j] > lo[j]) {
      tot <- tot + prim(hi[j], cc[j]) - prim(lo[j], cc[j])
    }
  }
  k * tot
}

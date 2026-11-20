#' Moment estimator of the extreme value index (Dekkers, Einmahl & de Haan
#' 1989)
#'
#' Estimates the extreme value index `gamma`, for any sign, from the first
#' two moments of the log-excesses of the `k` largest observations over the
#' `(k + 1)`-th largest, `X[n-k,n]`:
#' `M1 = mean(log(X[n-i,n]) - log(X[n-k,n]))` (the Hill estimator, equation
#' 1.3), `M2` the mean of the squares (1.5), and
#' `gamma = M1 + 1 - 1 / (2 (1 - M1^2 / M2))` (1.7). Logarithms are used, so
#' the `k + 1` largest observations must be positive: the estimator is scale
#' but not location invariant, and the paper suggests shifting the data if
#' needed (p. 1834).
#'
#' Under the second-order conditions of Theorem 3.1, `sqrt(k) (gamma_hat -
#' gamma)` is asymptotically normal with variance `v(gamma) = 1 + gamma^2`
#' for `gamma >= 0` and
#' `(1 - gamma)^2 (1 - 2 gamma) (4 - 8 (1 - 2 gamma) / (1 - 3 gamma) +
#' (5 - 11 gamma) (1 - 2 gamma) / ((1 - 3 gamma) (1 - 4 gamma)))` for
#' `gamma < 0` (Corollary 3.2, equation 3.10); `se` is `sqrt(v / k)` with
#' the estimate plugged in. (The matrix of equation 3.7 from which it is
#' derived prints the factor `(1 - 3 gamma) (1 - 4 gamma)` of `s22` in the
#' numerator instead of the denominator; equation 3.10 is correct.)
#'
#' For `gamma_hat < 0`, the end-point estimate is
#' `X[n-k,n] M1 (1 - 1 / gamma_hat) + X[n-k,n]` (equation 4.8). The paper
#' also gives an asymptotic variance for it (equation 5.3), but in our
#' simulations (generalized Pareto tails with `gamma = -0.5` and `-0.3`,
#' `n = 10^7`, `k = 400`) neither that variance nor a corrected delta-method
#' version matched the sampling variance, so no interval is returned. The
#' paper has no numerical example.
#'
#' @param x Numeric sample; the `k + 1` largest values must be positive.
#' @param k Number of log-excesses: whole numbers from 2 to `n - 1`, or
#'   `NULL` for all of them.
#'
#' @return An `evt_estimate` data frame with columns `k`, `gamma`, `se`,
#'   `M1`, `M2` and `endpoint` (`NA` unless `gamma < 0`).
#'
#' @references
#' Dekkers, A. L. M., Einmahl, J. H. J., & de Haan, L. (1989). A moment
#'   estimator for the index of an extreme-value distribution. The Annals of
#'   Statistics, 17(4), 1833-1855. doi:10.1214/aos/1176347397
#'
#' @examples
#' set.seed(1)
#' x <- 3 - stats::runif(2000)^0.5   # end-point 3, index -0.5
#' dekkers1989(x, k = c(50, 100, 200))
#'
#' @export
dekkers1989 <- function(x, k = NULL) {
  z <- evt_sorted(x, "Dekkers, Einmahl & de Haan (1989)", min_n = 3L)
  n <- length(z)
  k <- evt_k(k, 2L, n - 1L, "Dekkers, Einmahl & de Haan (1989)")
  if (any(z[seq_len(max(k) + 1L)] <= 0)) {
    stop("Dekkers, Einmahl & de Haan (1989) needs the k + 1 largest ",
         "observations to be positive (shift the data if needed).",
         call. = FALSE)
  }
  lz <- log(z[seq_len(max(k) + 1L)])
  rows <- lapply(k, function(kk) {
    d <- lz[seq_len(kk)] - lz[kk + 1L]
    m1 <- mean(d)
    m2 <- mean(d^2)
    g <- m1 + 1 - 0.5 / (1 - m1^2 / m2)
    end <- if (is.finite(g) && g < 0) {
      z[kk + 1L] * m1 * (1 - 1 / g) + z[kk + 1L]
    } else NA_real_
    data.frame(k = kk, gamma = g, se = sqrt(moment_variance(g) / kk),
               M1 = m1, M2 = m2, endpoint = end)
  })
  new_evt_estimate(do.call(rbind, rows), "Dekkers, Einmahl & de Haan (1989)")
}

#' Asymptotic variance of the moment estimator, equation 3.10
#' @noRd
moment_variance <- function(g) {
  if (!is.finite(g)) return(NA_real_)
  if (g >= 0) return(1 + g^2)
  (1 - g)^2 * (1 - 2 * g) * (4 - 8 * (1 - 2 * g) / (1 - 3 * g) +
                               (5 - 11 * g) * (1 - 2 * g) /
                               ((1 - 3 * g) * (1 - 4 * g)))
}

#' Neves & Pereira (2010) tests of a finite right end-point
#'
#' Tests `H0`: the distribution is in the Gumbel domain with an infinite
#' right end-point, against `H1`: it is in the Weibull or Gumbel domain
#' (`gamma <= 0`) with a finite end-point (Section 2). With the sample in
#' decreasing order and positive top observations:
#' * `T1 = mean_{i=1..k} (X[n-i,n] - X[n-k,n] - a) / (X[n,n] - X[n-k,n])`
#'   (equation 5), where `a = X[n-k,n] (M1 / 2) / (1 - M1^2 / M2)` is the
#'   moment estimator of the scale (equations 17-18), with `M1`, `M2` the
#'   first two moments of `log X[n-i+1,n] - log X[n-k,n]`, `i = 1, ..., k`;
#'   note that the sum in `T1` starts at the second largest observation while
#'   the moments include the largest, as printed;
#' * `T2 = (1/k) sum_{i=1}^{k-1} i (X[n-i+1,n] - X[n-i,n]) / X[n-k,n]`
#'   (equation 6).
#' The normalised statistics `T1* = sqrt(k) log(k) T1` and
#' `T2* = sqrt(k) (log(n/k) T2 - 1)` are asymptotically standard normal under
#' `H0` (Theorem 1 and Remark 3), and the paper rejects `H0` when `|T*|`
#' exceeds the normal quantile; `p_value` is the two-sided normal p-value.
#' `T2*` uses the sample size `n`. Neither statistic is location invariant
#' (`T1` uses logarithms through `a`, and `T2` divides by `X[n-k,n]`): a shift
#' of the data changes them, and the tests are calibrated for the
#' distributions as they are, for example the standard exponential under
#' `H0`, as in the paper's simulations.
#'
#' Two statements of the paper do not agree with our checks (see the
#' reproduction vignette): the limit in equation 9 is printed as
#' `B1 + B2`, while the construction of `T1` gives `B1 - B2`; and `T2*` is said
#' to shift to the right for distributions in the Gumbel domain with a finite
#' end-point, while for the paper's own negative Frechet example it shifts to
#' the left (both lead to rejection by the two-sided test).
#'
#' @param x Numeric sample; the `k + 1` largest values must be positive.
#' @param k Number of excesses: whole numbers from 2 to `n - 1`, or `NULL`
#'   for all of them.
#' @param statistic `"T1"` or `"T2"`.
#'
#' @return An `evt_test` data frame with columns `k`, `statistic` (`T1*` or
#'   `T2*`) and `p_value`.
#'
#' @references
#' Neves, C., & Pereira, A. (2010). Detecting finiteness in the right
#'   endpoint of light-tailed distributions. Statistics & Probability
#'   Letters, 80(5-6), 437-444. doi:10.1016/j.spl.2009.11.021
#'
#' @examples
#' set.seed(1)
#' neves_pereira2010(stats::rexp(1000) + 1, k = c(50, 100))
#' neves_pereira2010(2 - stats::runif(1000)^0.3, k = c(50, 100))
#'
#' @export
neves_pereira2010 <- function(x, k = NULL, statistic = c("T1", "T2")) {
  statistic <- match.arg(statistic)
  label <- "Neves & Pereira (2010)"
  z <- evt_sorted(x, label, min_n = 3L)
  n <- length(z)
  k <- evt_k(k, 2L, n - 1L, label)
  if (any(z[seq_len(max(k) + 1L)] <= 0)) {
    stop(label, " needs the k + 1 largest observations to be positive.",
         call. = FALSE)
  }
  rows <- lapply(k, function(kk) {
    thr <- z[kk + 1L]
    if (!(z[1L] > thr)) {
      stop("the k + 1 largest observations must not all be equal.",
           call. = FALSE)
    }
    if (statistic == "T1") {
      d <- log(z[seq_len(kk)]) - log(thr)
      m1 <- mean(d)
      m2 <- mean(d^2)
      a <- thr * (m1 / 2) / (1 - m1^2 / m2)
      t1 <- mean(z[seq_len(kk) + 1L] - thr - a) / (z[1L] - thr)
      s <- sqrt(kk) * log(kk) * t1
    } else {
      i <- seq_len(kk - 1L)
      t2 <- sum(i * (z[i] - z[i + 1L])) / (kk * thr)
      s <- sqrt(kk) * (log(n / kk) * t2 - 1)
    }
    data.frame(k = kk, statistic = s,
               p_value = 2 * stats::pnorm(-abs(s)))
  })
  new_evt_test(do.call(rbind, rows), paste0(label, ", ", statistic, "*"),
               "Gumbel domain with infinite end-point", "two.sided")
}

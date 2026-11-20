#' Neves, Picek & Fraga Alves (2006) ratio test of the Gumbel domain
#'
#' Tests `H0: gamma = 0` (the Gumbel domain) against `gamma != 0` with the
#' ratio of the largest excess to the mean excess of the `k` largest
#' observations over the `(k + 1)`-th largest:
#' `T = (X[n,n] - X[n-k,n]) / mean(X[n-i+1,n] - X[n-k,n])` (equation 3) and
#' `T* = T - log(k)` (equation 7). It is location and scale invariant. Under
#' `H0` and a second-order condition on `k` (Theorem 6), `T*` converges to
#' the standard Gumbel distribution; it tends to `-Inf` when `gamma < 0` and to
#' `+Inf` when `gamma > 0` (Theorem 9). The paper rejects with Gumbel
#' quantiles; `p_value` is the corresponding Gumbel p-value.
#'
#' For an exponential tail the excesses are, given the threshold, the order
#' statistics of `k` independent exponentials, so `T / k` is the ratio of
#' the largest to the sum of `k` exponentials, whose distribution is known
#' in closed form (Fisher 1929, the null distribution of his `g` statistic:
#' `P(max / sum > x) = sum_j (-1)^(j - 1) choose(k, j) (1 - j x)^(k - 1)`
#' over `j <= 1 / x`). `p_exact` gives the p-value under this exact null for
#' fixed `k` (Example 8 of the paper shows this representation), or `NA`
#' when the alternating sum has terms too large to be summed accurately
#' (for some statistics when `k` is in the thousands). The paper itself uses
#' the Gumbel limit.
#'
#' @param x Numeric sample.
#' @param k Number of excesses: whole numbers from 2 to `n - 1`, or `NULL`
#'   for all of them.
#' @param alternative `"two.sided"` (`gamma != 0`), `"positive"` (large `T*`)
#'   or `"negative"` (small `T*`).
#'
#' @return An `evt_test` data frame with columns `k`, `statistic` (`T*`),
#'   `p_value` (Gumbel limit) and `p_exact`.
#'
#' @references
#' Neves, C., Picek, J., & Fraga Alves, M. I. (2006). The contribution of the
#'   maximum to the sum of excesses for testing max-domains of attraction.
#'   Journal of Statistical Planning and Inference, 136(4), 1281-1301.
#'   doi:10.1016/j.jspi.2004.09.008
#'
#' Fisher, R. A. (1929). Tests of significance in harmonic analysis.
#'   Proceedings of the Royal Society of London A, 125(796), 54-59.
#'
#' @examples
#' set.seed(1)
#' neves2006(stats::rexp(1000), k = c(50, 100, 200))
#' neves2006(stats::runif(1000), k = c(50, 100, 200))
#'
#' @export
neves2006 <- function(x, k = NULL,
                      alternative = c("two.sided", "positive", "negative")) {
  alternative <- evt_alternative(alternative)
  label <- "Neves, Picek & Fraga Alves (2006)"
  z <- evt_sorted(x, label, min_n = 3L)
  n <- length(z)
  k <- evt_k(k, 2L, n - 1L, label)
  rows <- lapply(k, function(kk) {
    e <- evt_excess(z, kk)
    if (!(mean(e) > 0)) {
      stop("the k + 1 largest observations must not all be equal.",
           call. = FALSE)
    }
    tn <- e[1L] / mean(e)
    ts <- tn - log(kk)
    lam <- exp(-exp(-ts))                # Gumbel distribution function
    up <- fisher_g_upper(tn / kk, kk)    # P(T >= t) under the exact null
    data.frame(k = kk, statistic = ts,
               p_value = evt_pick(1 - lam, lam, alternative),
               p_exact = evt_pick(up, 1 - up, alternative))
  })
  new_evt_test(do.call(rbind, rows), label, "gamma = 0 (Gumbel domain)",
               alternative)
}

#' P(max / sum > x) for k independent exponentials (Fisher 1929)
#' @noRd
fisher_g_upper <- function(x, k) {
  if (x >= 1) return(0)
  if (x <= 1 / k) return(1)
  j <- seq_len(floor(1 / x))
  terms <- (-1)^(j - 1) * exp(lchoose(k, j) + (k - 1) * log1p(-j * x))
  # the alternating sum loses accuracy when its terms are very large
  if (!all(is.finite(terms)) || max(abs(terms)) > 1e10) return(NA_real_)
  min(1, max(0, sum(terms)))
}

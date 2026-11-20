#' Hasofer & Wang (1992) test of the Gumbel domain of attraction
#'
#' Tests `H0: gamma = 0` (the Gumbel domain) against `gamma != 0` with the
#' statistic of Hasofer & Wang (1992), computed on the `k + 1` largest
#' observations `X[1] >= ... >= X[k + 1]`:
#' `W = (k + 1) (mean(X) - X[k + 1])^2 / (k sum((X - mean(X))^2))`, the
#' squared correlation-type statistic of their Section 4. It is location and
#' scale invariant. Small values point to `gamma > 0` and large values to
#' `gamma < 0` (p. 173).
#'
#' **Convention for `k`.** Hasofer & Wang use the letter `k` for the number of
#' top order statistics, threshold included; here, as in the rest of the
#' package, `k` is the number of excesses over the threshold, so their
#' `k = 50` is `k = 49` here.
#'
#' **Null distribution.** Under `H0`, with `k` fixed and `n` large, the
#' excesses `X[i] - X[k + 1]` are distributed as the order statistics of `k`
#' independent exponentials, and `W = 1 / (k ((k + 1) G - 1))` with `G` the
#' Greenwood statistic `sum(E^2) / sum(E)^2` of `k` exponentials (Section 7).
#' `null = "simulated"` (default) uses this exact limiting null, simulated
#' with a fixed seed; `null = "normal"` uses the normal approximation of
#' Section 7, mean `1 / k` and variance `4 (k - 1) / (k^2 (k + 2) (k + 3))`
#' (from Moran 1947). The percentage points of their Table 8 were interpolated
#' from Stephens (1981); the simulated null reproduces the lower tail and the
#' upper 10% point, but the upper 5%, 2.5% and 1% points of the table are too
#' small for fewer than about 25 order statistics (an interpolation error in
#' the table).
#'
#' **Reproduction.** The simulations of Tables 1 and 2 (`n = 300`, 50 top
#' order statistics) are reproduced, as are the powers of Table 3 against
#' `gamma > 0`. Against `gamma < 0` the published powers are lower than those
#' obtained with the critical values of Table 2 (for `gamma = -1/2` at level
#' .01, .600 printed against about .69), although the distribution of `W`
#' agrees with Table 1; the paper does not state the critical values used
#' for those rows.
#'
#' @param x Numeric sample.
#' @param k Number of excesses (the statistic uses the `k + 1` largest
#'   values): whole numbers from 2 to `n - 1`, or `NULL` for all of them.
#' @param alternative `"two.sided"` (`gamma != 0`), `"positive"`
#'   (`gamma > 0`, small `W`) or `"negative"` (`gamma < 0`, large `W`).
#' @param null `"simulated"` or `"normal"`, see above.
#' @param nsim Number of simulations of the null distribution.
#'
#' @return An `evt_test` data frame with columns `k`, `statistic` (`W`) and
#'   `p_value`.
#'
#' @references
#' Hasofer, A. M., & Wang, Z. (1992). A test for extreme value domain of
#'   attraction. Journal of the American Statistical Association, 87(417),
#'   171-177. doi:10.1080/01621459.1992.10475189
#'
#' @seealso [neves_fraga_alves2007()] for the version with an intermediate
#'   number of excesses.
#'
#' @examples
#' set.seed(1)
#' hasofer_wang1992(stats::rexp(300), k = c(19, 49))    # Gumbel domain
#' hasofer_wang1992(stats::runif(300), k = 49)          # gamma = -1
#'
#' @export
hasofer_wang1992 <- function(x, k = NULL,
                             alternative = c("two.sided", "positive",
                                             "negative"),
                             null = c("simulated", "normal"), nsim = 1e5) {
  alternative <- evt_alternative(alternative)
  null <- match.arg(null)
  check_nsim(nsim)
  z <- evt_sorted(x, "Hasofer & Wang (1992)", min_n = 3L)
  n <- length(z)
  k <- evt_k(k, 2L, n - 1L, "Hasofer & Wang (1992)")
  rows <- lapply(k, function(kk) {
    w <- hw_stat(z[seq_len(kk + 1L)])
    if (null == "simulated") {
      nul <- hw_null(kk, nsim)
      p_pos <- mean(nul <= w)
      p_neg <- mean(nul >= w)
    } else {
      s <- sqrt(4 * (kk - 1) / (kk^2 * (kk + 2) * (kk + 3)))
      p_pos <- stats::pnorm((w - 1 / kk) / s)
      p_neg <- 1 - p_pos
    }
    data.frame(k = kk, statistic = w,
               p_value = evt_pick(p_pos, p_neg, alternative))
  })
  new_evt_test(do.call(rbind, rows), "Hasofer & Wang (1992)",
               "gamma = 0 (Gumbel domain)", alternative)
}

#' W of Hasofer & Wang on the top order statistics (decreasing)
#' @noRd
hw_stat <- function(x) {
  m <- length(x)
  xbar <- mean(x)
  ss <- sum((x - xbar)^2)
  if (!(ss > 0)) {
    stop("the k + 1 largest observations must not all be equal.",
         call. = FALSE)
  }
  m * (xbar - x[m])^2 / ((m - 1) * ss)
}

#' Limiting null distribution of W for k excesses (k exponentials)
#' @noRd
hw_null <- function(k, nsim) {
  sim_null(paste("hw", k), nsim, function(ns) {
    g <- greenwood_sim(k, ns)
    1 / (k * ((k + 1) * g - 1))
  })
}

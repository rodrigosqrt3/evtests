#' Neves & Fraga Alves (2007) tests of the Gumbel domain with an intermediate
#' number of excesses
#'
#' Tests `H0: gamma = 0` (the Gumbel domain) against `gamma != 0` with the
#' standardised Greenwood-type and Hasofer-Wang statistics of Neves & Fraga
#' Alves (2007), for `k = k(n)` intermediate (`k -> Inf`, `k / n -> 0`). With
#' the excesses `Z[i] = X[n-i+1,n] - X[n-k,n]`, `i = 1, ..., k`, and their
#' first two moments `M1`, `M2`:
#' * `R = M2 / M1^2` (equation 1.3) and `R* = sqrt(k / 4) (R - 2)` (3.1);
#' * `W* = sqrt(k / 4) (k W - 1)` (3.2), with `W` the Hasofer-Wang statistic.
#' Both are asymptotically standard normal under `H0` when
#' `sqrt(k) A(n/k) -> 0` (Theorem 3.1 and Remark 3.1). Large `R*` points to
#' `gamma > 0` and small `R*` to `gamma < 0`; for `W*` it is the reverse
#' (p. 304). p-values are from the normal distribution.
#'
#' **Which W.** The paper prints `W = 1 / (k (R - 1))` (equations 1.2 and
#' 2.3), but its Table 1 is reproduced only with the original statistic of
#' Hasofer & Wang (1992) on the `k + 1` largest observations,
#' `W = 1 / ((k + 1) R - k)` (see [hasofer_wang1992()]); with the printed
#' formula the size at `gamma = 0`, `k = 50` is about 0.076 instead of the
#' published 0.051. This function uses the original statistic. The two
#' agree asymptotically.
#'
#' @param x Numeric sample.
#' @param k Number of excesses: whole numbers from 2 to `n - 1`, or `NULL`
#'   for all of them.
#' @param statistic `"R"` for `R*` (the paper's "Gt" test) or `"W"` for `W*`
#'   (its "HW" test).
#' @param alternative `"two.sided"` (`gamma != 0`), `"positive"` or
#'   `"negative"` (the sign of `gamma` under the alternative).
#'
#' @return An `evt_test` data frame with columns `k`, `statistic` (`R*` or
#'   `W*`) and `p_value`.
#'
#' @references
#' Neves, C., & Fraga Alves, M. I. (2007). Semi-parametric approach to the
#'   Hasofer-Wang and Greenwood statistics in extremes. TEST, 16(2), 297-313.
#'   doi:10.1007/s11749-006-0010-1
#'
#' @examples
#' set.seed(1)
#' x <- stats::rexp(1000)
#' neves_fraga_alves2007(x, k = c(50, 100, 200))
#' neves_fraga_alves2007(x, k = c(50, 100, 200), statistic = "W")
#'
#' @export
neves_fraga_alves2007 <- function(x, k = NULL, statistic = c("R", "W"),
                                  alternative = c("two.sided", "positive",
                                                  "negative")) {
  statistic <- match.arg(statistic)
  alternative <- evt_alternative(alternative)
  label <- "Neves & Fraga Alves (2007)"
  z <- evt_sorted(x, label, min_n = 3L)
  n <- length(z)
  k <- evt_k(k, 2L, n - 1L, label)
  rows <- lapply(k, function(kk) {
    e <- evt_excess(z, kk)
    m1 <- mean(e)
    if (!(m1 > 0)) {
      stop("the k + 1 largest observations must not all be equal.",
           call. = FALSE)
    }
    r <- mean(e^2) / m1^2
    if (statistic == "R") {
      s <- sqrt(kk / 4) * (r - 2)
      p_pos <- 1 - stats::pnorm(s)
    } else {
      w <- 1 / ((kk + 1) * r - kk)
      s <- sqrt(kk / 4) * (kk * w - 1)
      p_pos <- stats::pnorm(s)
    }
    data.frame(k = kk, statistic = s,
               p_value = evt_pick(p_pos, 1 - p_pos, alternative))
  })
  new_evt_test(do.call(rbind, rows),
               paste0(label, ", ", statistic, "*"),
               "gamma = 0 (Gumbel domain)", alternative)
}

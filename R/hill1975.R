#' Hill (1975) estimator of a Pareto-type tail
#'
#' Conditional maximum likelihood estimators of Hill (1975) for an upper tail
#' of the form `1 - F(y) = C y^(-alpha)` beyond some point, computed from the
#' `k + 1` largest observations `y[1] >= ... >= y[k + 1]` (Hill's `r` is `k`
#' here). With the spacings `V[i] = log(y[i] / y[i + 1])`,
#' * `alpha0 = k / sum(i * V[i]) = k / (sum(log(y[1:k])) - k log(y[k + 1]))`
#'   (p. 1166 and p. 1172), whose inverse `gamma = 1 / alpha0` is the usual
#'   Hill estimator of the extreme value index;
#' * `alpha1 = (k + 1) / k * alpha0` and `C1 = y[k + 1]^alpha1 (k + 1) / n`
#'   (equations 2.10 and 3.1), the estimates from the likelihood that also
#'   uses the position of `y[k + 1]`.
#'
#' Given that the `k + 1` largest observations lie in the Pareto part of the
#' tail, `alpha * i * V[i]` are independent standard exponential variables, so
#' `alpha0` is distributed as `k alpha / G` with `G` a Gamma(k, 1) variable,
#' with mean `k alpha / (k - 1)` and standard deviation `alpha` times
#' `sd_rel = k / ((k - 1) sqrt(k - 2))` (equation 2.6). The paper gives no
#' interval formula. To judge how many observations follow the Pareto tail,
#' Hill proposes to compare the `i V[i]` with an exponential sample, and the
#' discrepancy statistic `H = alpha0^2 sum((i V[i] - 1 / alpha0)^2)`, which
#' is about `k - 1` on average with variance about `8 (k - 1)` (p. 1172).
#'
#' @param x Numeric sample; the `k + 1` largest values must be positive.
#' @param k Number of spacings (the estimators use the `k + 1` largest
#'   values): whole numbers from 1 to `n - 1`, or `NULL` for all of them.
#'
#' @return An `evt_estimate` data frame with one row per `k` and columns `k`,
#'   `gamma` (`1 / alpha0`), `alpha0`, `alpha1`, `C1`, `sd_rel` (`NA` for
#'   `k <= 2`) and `H` (0 for `k = 1`), plus the attribute `iV`, the values
#'   `i V[i]` for `i = 1, ..., max(k)`.
#'
#' @references
#' Hill, B. M. (1975). A simple general approach to inference about the tail
#'   of a distribution. The Annals of Statistics, 3(5), 1163-1174.
#'   doi:10.1214/aos/1176343247
#'
#' @seealso [us_cities_1940], the example of the paper.
#'
#' @examples
#' # Table 1 of Hill (1975): the 30 largest cities of the United States in 1940
#' h <- hill1975(us_cities_1940$population, k = 1:29)
#' h[c(1, 10, 29), c("k", "alpha0", "alpha1", "H")]
#' plot(h, "alpha0")
#'
#' @export
hill1975 <- function(x, k = NULL) {
  y <- evt_sorted(x, "Hill (1975)")
  n <- length(y)
  k <- evt_k(k, 1L, n - 1L, "Hill (1975)")
  kmax <- max(k)
  if (any(y[seq_len(kmax + 1L)] <= 0)) {
    stop("Hill (1975) needs the k + 1 largest observations to be positive.",
         call. = FALSE)
  }
  ly <- log(y[seq_len(kmax + 1L)])
  iv <- seq_len(kmax) * (ly[seq_len(kmax)] - ly[seq_len(kmax) + 1L])
  rows <- lapply(k, function(r) {
    s <- sum(ly[seq_len(r)]) - r * ly[r + 1L]
    a0 <- r / s
    a1 <- (r + 1) / r * a0
    c1 <- exp(a1 * ly[r + 1L]) * (r + 1) / n
    sd_rel <- if (r > 2L) r / ((r - 1) * sqrt(r - 2)) else NA_real_
    h <- a0^2 * sum((iv[seq_len(r)] - 1 / a0)^2)
    data.frame(k = r, gamma = 1 / a0, alpha0 = a0, alpha1 = a1, C1 = c1,
               sd_rel = sd_rel, H = h)
  })
  out <- new_evt_estimate(do.call(rbind, rows), "Hill (1975)")
  attr(out, "iV") <- iv
  out
}

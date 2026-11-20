#' Pickands (1975) estimator of the extreme value index
#'
#' Estimates the shape `c` (the extreme value index) and the scale `a` of the
#' generalized Pareto approximation to the excesses over a high level, from
#' three upper order statistics. With the sample in decreasing order
#' `Z[1] >= Z[2] >= ...` and an integer `M` with `4M <= n` (`k` here),
#' * `c = log((Z[M] - Z[2M]) / (Z[2M] - Z[4M])) / log(2)` (equation 2.3),
#' * `a = (Z[2M] - Z[4M]) / integral_0^log(2) exp(c s) ds`, that is
#'   `(Z[2M] - Z[4M]) c / (2^c - 1)`, or `(Z[2M] - Z[4M]) / log(2)` when
#'   `c = 0` (equation 2.4; equation 4.3 prints `exp(-c s)`, a misprint, since
#'   equations 2.4 and 3.3 and p. 131 have `exp(c s)`).
#' These follow from equating the 1/2 and 1/4 percentiles of the generalized
#' Pareto distribution of the excesses over `Z[4M]` to the empirical ones
#' (equations 3.1-3.3). Both estimates are location and scale equivariant.
#'
#' Pickands proposes to choose `M` by minimising the largest distance between
#' the empirical and the fitted distribution of the excesses, but is not
#' consistent about which of several minimisers to take (p. 122 and p. 130),
#' so this function leaves `k` to the user; the paper gives no variance, no
#' interval and no numerical example.
#'
#' @param x Numeric sample.
#' @param k The integer `M`: whole numbers with `4 k <= n`, or `NULL` for all
#'   of them.
#'
#' @return An `evt_estimate` data frame with columns `k` (`M`), `gamma` (`c`)
#'   and `scale` (`a`). `gamma` is not finite when two of the order statistics
#'   used are tied.
#'
#' @references
#' Pickands, J. (1975). Statistical inference using extreme order statistics.
#'   The Annals of Statistics, 3(1), 119-131. doi:10.1214/aos/1176343003
#'
#' @examples
#' set.seed(1)
#' x <- rexp(400)          # exponential tail: index 0
#' pickands1975(x, k = c(10, 25, 50, 100))
#'
#' @export
pickands1975 <- function(x, k = NULL) {
  z <- evt_sorted(x, "Pickands (1975)", min_n = 4L)
  n <- length(z)
  k <- evt_k(k, 1L, n %/% 4L, "Pickands (1975)")
  d1 <- z[k] - z[2L * k]
  d2 <- z[2L * k] - z[4L * k]
  g <- log(d1 / d2) / log(2)
  integral <- ifelse(abs(g) < 1e-12, log(2), (2^g - 1) / g)
  new_evt_estimate(data.frame(k = k, gamma = g, scale = d2 / integral),
                   "Pickands (1975)")
}

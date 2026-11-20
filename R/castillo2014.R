#' del Castillo, Daoudi & Lockhart (2014) tests of an exponential tail based
#' on the residual coefficient of variation
#'
#' Tests whether the excesses over high thresholds are exponential (the
#' residual coefficient of variation equals 1, as in the Gumbel domain with an
#' exponential-type tail) against heavier or lighter tails. For a threshold
#' `t`, the residual coefficient of variation `cv(t)` is the standard
#' deviation (divisor `n(t)`) of the observations above `t` over the mean of
#' their excesses `X - t` (equation 1); under exponentiality
#' `sqrt(n(t)) (cv(t) - 1)` is asymptotically standard normal (equation 7).
#' The test statistic combines the thresholds `q_k` that leave the
#' `n_k = n / 2^k` largest observations above them, `k = 0, ..., m`
#' (`q_0 = 0`, the whole sample):
#' `T_m = sum_{k=0}^m n_k (cv(q_k) - 1)^2` (equations 9-10), large values
#' rejecting. The data must be non-negative.
#'
#' `p_value` is computed, as in the paper's applications, from simulated
#' exponential samples of the same size (`nsim` samples, fixed seed).
#' `p_approx` uses the three-moment approximation `a + b chi^2_nu` to the
#' asymptotic distribution of Theorem 2, the quadratic form
#' `sum lambda_i Z_i^2` with `lambda` the eigenvalues of the matrix with
#' entries `2^(-|i - j| / 2)` (Section 4.2): `b = sum lambda^3 /
#' sum lambda^2`, `nu = (sum lambda^2)^3 / (sum lambda^3)^2`,
#' `a = sum lambda - b nu`, which reproduces the "Approximate" row of the
#' paper's Table 1.
#'
#' The paper does not say how `n / 2^k` is rounded when it is not an integer;
#' this function rounds down. With this choice the finite-sample 95% points of
#' the paper's Table 1 are reproduced within about 0.4 (0.1 for `n >= 500`),
#' but the 99% points of `T_1` to `T_4` differ by up to 1.8 in both
#' directions; the printed 99% points are not monotone in `n` (see the
#' reproduction vignette).
#'
#' @param x Non-negative numeric sample.
#' @param m Number of halvings of the sample, a whole number of at least 0
#'   with `n / 2^m >= 2`.
#' @param nsim Number of simulations of the null distribution.
#'
#' @return An `evt_test` data frame with one row and columns `m`,
#'   `statistic` (`T_m`), `p_value` and `p_approx`.
#'
#' @references
#' del Castillo, J., Daoudi, J., & Lockhart, R. (2014). Methods to
#'   distinguish between polynomial and exponential tails. Scandinavian
#'   Journal of Statistics, 41(2), 382-393. doi:10.1111/sjos.12037
#'
#' @seealso [castillo2014_cv()] for the CV-plot.
#'
#' @examples
#' set.seed(1)
#' castillo2014(stats::rexp(500), m = 4, nsim = 2000)
#' castillo2014(stats::rexp(500)^1.5, m = 4, nsim = 2000)
#'
#' @export
castillo2014 <- function(x, m = 4L, nsim = 1e4) {
  check_nsim(nsim)
  label <- "del Castillo, Daoudi & Lockhart (2014)"
  z <- evt_sorted(x, label, min_n = 4L)
  if (any(z < 0)) stop(label, " needs non-negative data.", call. = FALSE)
  n <- length(z)
  if (!is.numeric(m) || length(m) != 1L || m < 0 || m != round(m) ||
      n %/% 2^m < 2) {
    stop("`m` must be a whole number of at least 0 with n / 2^m >= 2.",
         call. = FALSE)
  }
  s <- castillo_stat(z, m)
  nul <- sim_null(paste("castillo", n, m), nsim, function(ns) {
    vapply(seq_len(ns), function(i) {
      castillo_stat(sort(stats::rexp(n), decreasing = TRUE), m)
    }, numeric(1))
  })
  ap <- castillo_approx(m)
  p_ap <- 1 - stats::pchisq(max(0, (s - ap[["a"]]) / ap[["b"]]), ap[["nu"]])
  new_evt_test(data.frame(m = m, statistic = s, p_value = mean(nul >= s),
                          p_approx = p_ap),
               paste0(label, ", T", m), "exponential excesses (cv = 1)",
               "two.sided")
}

#' Residual coefficient of variation of the j largest observations over the
#' next one (z decreasing; j = n uses the whole sample and threshold 0)
#' @noRd
residual_cv <- function(z, j) {
  top <- z[seq_len(j)]
  thr <- if (j >= length(z)) 0 else z[j + 1L]
  sqrt(mean(top^2) - mean(top)^2) / mean(top - thr)
}

#' T_m of equation 10 with n_k = floor(n / 2^k)
#' @noRd
castillo_stat <- function(z, m) {
  n <- length(z)
  nk <- n %/% 2^(0:m)
  sum(vapply(nk, function(j) j * (residual_cv(z, j) - 1)^2, numeric(1)))
}

#' Three-moment approximation a + b chi^2_nu to the limit of T_m
#' @noRd
castillo_approx <- function(m) {
  idx <- 0:m
  lam <- eigen(2^(-abs(outer(idx, idx, "-")) / 2), symmetric = TRUE,
               only.values = TRUE)$values
  b <- sum(lam^3) / sum(lam^2)
  nu <- sum(lam^2)^3 / sum(lam^3)^2
  c(a = sum(lam) - b * nu, b = b, nu = nu)
}

#' CV-plot of del Castillo, Daoudi & Lockhart (2014)
#'
#' The residual coefficient of variation of the excesses of the `k` largest
#' observations over the `(k + 1)`-th largest, as a function of `k`, with the
#' pointwise limits `1 +/- z / sqrt(k)` under exponentiality (equation 7 and
#' Section 3; the paper draws 90% limits and uses the upper one, since
#' Pareto-type tails give `cv > 1`). `p_value` is the one-sided upper normal
#' p-value of `sqrt(k) (cv - 1)`.
#'
#' @param x Non-negative numeric sample.
#' @param k Number of excesses, whole numbers from 2 to `n - 1`, or `NULL`
#'   for all of them.
#' @param level Level of the pointwise limits.
#'
#' @return An `evt_estimate` data frame with columns `k`, `cv`, `lower`,
#'   `upper` and `p_value`; `plot()` draws the CV-plot.
#'
#' @references
#' del Castillo, J., Daoudi, J., & Lockhart, R. (2014). Methods to
#'   distinguish between polynomial and exponential tails. Scandinavian
#'   Journal of Statistics, 41(2), 382-393. doi:10.1111/sjos.12037
#'
#' @examples
#' set.seed(1)
#' cv <- castillo2014_cv(stats::rexp(1000), k = 10:900)
#' plot(cv, "cv")
#'
#' @export
castillo2014_cv <- function(x, k = NULL, level = 0.9) {
  check_alpha(level)
  label <- "del Castillo, Daoudi & Lockhart (2014)"
  z <- evt_sorted(x, label, min_n = 3L)
  if (any(z < 0)) stop(label, " needs non-negative data.", call. = FALSE)
  n <- length(z)
  k <- evt_k(k, 2L, n - 1L, label)
  cv <- vapply(k, function(j) residual_cv(z, j), numeric(1))
  q <- stats::qnorm(1 - (1 - level) / 2)
  new_evt_estimate(data.frame(k = k, cv = cv, lower = 1 - q / sqrt(k),
                              upper = 1 + q / sqrt(k),
                              p_value = 1 - stats::pnorm(sqrt(k) * (cv - 1))),
                   paste0(label, ", CV-plot"), param = "cv")
}

#' Fraga Alves & Neves (2014) estimator of a finite right end-point
#'
#' Estimates the finite right end-point `x_F` of a distribution in the Gumbel
#' or Weibull domain of attraction (`gamma <= 0`) without estimating the
#' extreme value index. With the sample in decreasing order
#' `X[1] >= X[2] >= ...`,
#' `x_F = X[1] + X[k + 1] - (1 / log(2)) sum_{i=0}^{k-1}
#' log(1 + 1 / (k + i)) X[k + i + 1]`, equivalently
#' `X[1] + sum_{i=0}^{k-1} a[i] (X[k + 1] - X[k + i + 1])` with positive
#' weights `a[i] = log((k + i + 1) / (k + i)) / log(2)` summing to one
#' (Fraga Alves & Neves 2014; Fraga Alves, Neves & Rosario 2016, equations 6
#' and 7). It uses the `2k` largest observations and is never below the
#' sample maximum. It is consistent for `gamma <= 0` and diverges for
#' `gamma > 0` (Fraga Alves et al. 2016, Proposition 1).
#'
#' Fraga Alves, Neves & Rosario (2016) also give, for `-1/2 < gamma < 0`, an
#' approximate one-sided upper confidence bound (their equation 14):
#' `x_F - a (h(gamma) + k^gamma q)`, with
#' `h(gamma) = ((2^(-gamma) - 1) / (gamma log(2)) + 1) / gamma` (equation
#' 11) and `q = (-log(alpha))^(-gamma) / gamma`, where `gamma` and the scale
#' `a` are the estimators of their equations 19-21 from the `k` excesses
#' `e[i] = X[i] - X[k + 1]`: `N1 = mean(e)`, `N2 = mean(e^2)`,
#' `gamma = 1 - 1 / (2 (1 - N1^2 / N2))` and `a = N1 (1 - gamma)`. The same
#' `k` is used throughout, as in the paper's application. When the estimated
#' `gamma` is outside `(-1/2, 0)` the bound is `NA` with a warning. The
#' convergence to the limit is slow when `gamma` is close to 0, as the paper
#' notes; see the reproduction vignette for simulated coverages.
#'
#' **Interval in the Gumbel domain.** For a distribution in the Gumbel domain
#' with a finite end-point, Theorem 3 of Fraga Alves & Neves (2014) gives
#' `(x_F_hat - x_F) / a_hat -> Lambda - log(2) / 2 - lambda / log(2)`, with
#' `Lambda` standard Gumbel, `a_hat` the mean excess of the `k` largest
#' observations over `X[k + 1]` (which the paper calls the maximum likelihood
#' estimator of the generalized Pareto scale; it is that of the exponential
#' distribution) and `lambda` a bias term. With `interval = TRUE` the
#' equal-tailed interval
#' `[x_F_hat - a_hat (g(1 - alpha/2) - log(2)/2),
#' x_F_hat - a_hat (g(alpha/2) - log(2)/2)]`, `g` the Gumbel quantile
#' function, is returned, which assumes `lambda = 0`. In simulations of the
#' model of the paper's Section 5 (negative Frechet with `beta = 0.5`,
#' `n = 1000`, `k = 48`) the nominal 95% interval covered the end-point in all
#' of 500 samples: it is conservative at this sample size.
#'
#' @param x Numeric sample.
#' @param k Whole numbers from 1 to `n / 2`, or `NULL` for all of them.
#' @param bound If `TRUE`, add the upper bound of Fraga Alves et al. (2016);
#'   `k` must then be at least 2.
#' @param alpha Significance level of the bound and of the interval, in
#'   (0, 1).
#' @param interval If `TRUE`, add the interval of Theorem 3 of Fraga Alves &
#'   Neves (2014) for the Gumbel domain.
#'
#' @return An `evt_estimate` data frame with columns `k` and `endpoint`;
#'   with `bound = TRUE` also `gamma` and `scale` (equations 19-20) and
#'   `upper`; with `interval = TRUE` also `mean_excess`, `lower_gumbel` and
#'   `upper_gumbel`.
#'
#' @references
#' Fraga Alves, I., & Neves, C. (2014). Estimation of the finite right
#'   endpoint in the Gumbel domain. Statistica Sinica, 24(4), 1811-1835.
#'   doi:10.5705/ss.2013.183
#'
#' Fraga Alves, I., Neves, C., & Rosario, P. (2016). A general estimator for
#'   the right endpoint with an application to supercentenarian women's
#'   records. arXiv:1412.3972 (version 3).
#'
#' @examples
#' set.seed(1)
#' x <- 1 - stats::runif(1000)^0.4    # end-point 1, gamma = -0.4
#' fraga_alves_neves2014(x, k = c(25, 50, 100), bound = TRUE)
#'
#' @export
fraga_alves_neves2014 <- function(x, k = NULL, bound = FALSE, alpha = 0.05,
                                  interval = FALSE) {
  label <- "Fraga Alves & Neves (2014)"
  if (isTRUE(bound) || isTRUE(interval)) check_alpha(alpha)
  z <- evt_sorted(x, label)
  n <- length(z)
  k <- evt_k(k, if (isTRUE(bound)) 2L else 1L, n %/% 2L, label)
  rows <- lapply(k, function(kk) {
    out <- data.frame(k = kk, endpoint = fan_endpoint(z, kk))
    if (isTRUE(bound)) {
      e <- evt_excess(z, kk)
      n1 <- mean(e)
      n2 <- mean(e^2)
      g <- 1 - 0.5 / (1 - n1^2 / n2)
      a <- n1 * (1 - g)
      up <- NA_real_
      if (is.finite(g) && g > -0.5 && g < 0) {
        h <- ((2^(-g) - 1) / (g * log(2)) + 1) / g
        q <- (-log(alpha))^(-g) / g
        up <- out$endpoint - a * (h + kk^g * q)
      }
      out <- cbind(out, gamma = g, scale = a, upper = up)
    }
    if (isTRUE(interval)) {
      ah <- mean(evt_excess(z, kk))
      gq <- -log(-log(c(1 - alpha / 2, alpha / 2))) - log(2) / 2
      out <- cbind(out, mean_excess = ah,
                   lower_gumbel = out$endpoint - ah * gq[1L],
                   upper_gumbel = out$endpoint - ah * gq[2L])
    }
    out
  })
  out <- do.call(rbind, rows)
  if (isTRUE(bound) && anyNA(out$upper)) {
    warning("the estimated extreme value index is not in (-1/2, 0) for k = ",
            paste(out$k[is.na(out$upper)], collapse = ", "),
            ", where the bound of Fraga Alves, Neves & Rosario (2016) ",
            "holds; `upper` is NA there.", call. = FALSE)
  }
  new_evt_estimate(out, label, param = "endpoint")
}

#' Equation 6 of Fraga Alves et al. (2016), z decreasing, 2k values used
#' @noRd
fan_endpoint <- function(z, k) {
  i <- 0:(k - 1L)
  a <- (log(k + i + 1) - log(k + i)) / log(2)
  z[1L] + sum(a * (z[k + 1L] - z[k + i + 1L]))
}

#' Fraga Alves, Neves & Rosario (2016) test of the Gumbel domain
#'
#' Tests `H0: gamma = 0` (the Gumbel domain) against `gamma != 0` with the
#' statistic of Fraga Alves, Neves & Rosario (2016, equations 16-17), built on
#' the end-point estimator of [fraga_alves_neves2014()]:
#' `G = (x_F - X[k + 1]) / (X[k + 1] - X[2k + 1])` and
#' `G* = log(2) G - (log(k) + log(2) / 2)`, with the sample in decreasing
#' order (`2k + 1` largest observations used). Under `H0`, `G*` converges to
#' the standard Gumbel distribution; it tends to `+Inf` when `gamma > 0` and
#' to `-Inf` when `gamma < 0` (Theorem 6). `p_value` uses this limit, as the
#' paper does.
#'
#' For an exponential tail, given `X[2k + 1]` the `2k` larger observations
#' are distributed as `X[2k + 1]` plus the order statistics of `2k`
#' independent exponentials, so `G*` has a null distribution free of
#' parameters for fixed `k`. `p_sim` uses that distribution, simulated with
#' a fixed seed.
#'
#' @param x Numeric sample.
#' @param k Whole numbers from 1 to `(n - 1) / 2`, or `NULL` for all of them.
#' @param alternative `"two.sided"` (`gamma != 0`), `"positive"` (large `G*`)
#'   or `"negative"` (small `G*`).
#' @param nsim Number of simulations of the exponential-tail null.
#'
#' @return An `evt_test` data frame with columns `k`, `statistic` (`G*`),
#'   `p_value` (Gumbel limit) and `p_sim`.
#'
#' @references
#' Fraga Alves, I., Neves, C., & Rosario, P. (2016). A general estimator for
#'   the right endpoint with an application to supercentenarian women's
#'   records. arXiv:1412.3972 (version 3).
#'
#' @examples
#' set.seed(1)
#' fraga_alves2016(stats::rexp(1000), k = c(25, 50, 100))
#'
#' @export
fraga_alves2016 <- function(x, k = NULL,
                            alternative = c("two.sided", "positive",
                                            "negative"),
                            nsim = 1e4) {
  alternative <- evt_alternative(alternative)
  check_nsim(nsim)
  label <- "Fraga Alves, Neves & Rosario (2016)"
  z <- evt_sorted(x, label, min_n = 3L)
  n <- length(z)
  k <- evt_k(k, 1L, (n - 1L) %/% 2L, label)
  rows <- lapply(k, function(kk) {
    if (!(z[kk + 1L] > z[2L * kk + 1L])) {
      stop("the (k + 1)-th and the (2k + 1)-th largest observations must ",
           "differ.", call. = FALSE)
    }
    g <- fanr_stat(z, kk)
    lam <- exp(-exp(-g))
    nul <- sim_null(paste("fanr", kk), nsim, function(ns) {
      vapply(seq_len(ns), function(i) {
        fanr_stat(c(sort(stats::rexp(2L * kk), decreasing = TRUE), 0), kk)
      }, numeric(1))
    })
    data.frame(k = kk, statistic = g,
               p_value = evt_pick(1 - lam, lam, alternative),
               p_sim = evt_pick(mean(nul >= g), mean(nul <= g), alternative))
  })
  new_evt_test(do.call(rbind, rows), label, "gamma = 0 (Gumbel domain)",
               alternative)
}

#' G* of equation 17, z decreasing with at least 2k + 1 values
#' @noRd
fanr_stat <- function(z, k) {
  g <- (fan_endpoint(z, k) - z[k + 1L]) / (z[k + 1L] - z[2L * k + 1L])
  log(2) * g - (log(k) + log(2) / 2)
}

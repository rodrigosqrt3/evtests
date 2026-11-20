#' Maximum likelihood estimator of the extreme value index (Drees, Ferreira &
#' de Haan 2004)
#'
#' Fits the generalized Pareto distribution by maximum likelihood to the `k`
#' excesses `y[i] = X[n-i+1,n] - X[n-k,n]`, `i = 1, ..., k`, of the `k`
#' largest observations over the `(k + 1)`-th largest, and returns the
#' estimates of the extreme value index `gamma` and of the scale `sigma`.
#' The likelihood is maximised with [grimshaw1993()], to which the paper
#' defers the computation.
#'
#' Drees, Ferreira & de Haan (2004) prove that, for `gamma > -1/2` and an
#' intermediate sequence `k`, `sqrt(k) (gamma_hat - gamma, sigma_hat /
#' a(n/k) - 1)` is asymptotically normal with covariance matrix
#' `Sigma = [(1 + gamma)^2, -(1 + gamma); -(1 + gamma), 2 + 2 gamma +
#' gamma^2]` (Corollary 2.1), plus a bias term when `sqrt(k) A(n/k)` does not
#' vanish. The returned standard errors are `sqrt(diag(Sigma) / k)` with the
#' estimate plugged in, the scale one relative to `sigma`. The paper defines
#' the estimator through a maximum over `gamma > -1/2` and does not say what to
#' report otherwise; this function returns the unrestricted maximum
#' likelihood estimate and sets the standard errors to `NA` when
#' `gamma_hat <= -1/2`, where the result does not apply. The paper has no
#' numerical example.
#'
#' @param x Numeric sample.
#' @param k Number of excesses: whole numbers from 2 to `n - 1`, or `NULL`
#'   for all of them.
#'
#' @return An `evt_estimate` data frame with columns `k`, `gamma`, `scale`,
#'   `se_gamma`, `se_scale_rel` and `exists` (`FALSE` when no maximum
#'   likelihood estimate exists, with the estimates `NA`).
#'
#' @references
#' Drees, H., Ferreira, A., & de Haan, L. (2004). On maximum likelihood
#'   estimation of the extreme value index. The Annals of Applied
#'   Probability, 14(3), 1179-1201. doi:10.1214/105051604000000279
#'
#' Grimshaw, S. D. (1993). Computing maximum likelihood estimates for the
#'   generalized Pareto distribution. Technometrics, 35(2), 185-191.
#'
#' @examples
#' set.seed(1)
#' x <- stats::rexp(500)
#' drees2004(x, k = c(25, 50, 100))
#'
#' @export
drees2004 <- function(x, k = NULL) {
  z <- evt_sorted(x, "Drees, Ferreira & de Haan (2004)", min_n = 3L)
  n <- length(z)
  k <- evt_k(k, 2L, n - 1L, "Drees, Ferreira & de Haan (2004)")
  rows <- lapply(k, function(kk) {
    y <- z[seq_len(kk)] - z[kk + 1L]
    fit <- if (any(y > 0)) grimshaw1993(y) else list(exists = FALSE)
    if (!fit$exists) {
      return(data.frame(k = kk, gamma = NA_real_, scale = NA_real_,
                        se_gamma = NA_real_, se_scale_rel = NA_real_,
                        exists = FALSE))
    }
    g <- fit$gamma
    ok <- g > -0.5
    data.frame(k = kk, gamma = g, scale = fit$scale,
               se_gamma = if (ok) (1 + g) / sqrt(kk) else NA_real_,
               se_scale_rel = if (ok) sqrt((2 + 2 * g + g^2) / kk) else
                 NA_real_,
               exists = TRUE)
  })
  new_evt_estimate(do.call(rbind, rows), "Drees, Ferreira & de Haan (2004)")
}

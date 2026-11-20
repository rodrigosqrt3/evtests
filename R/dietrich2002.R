#' Tests of the extreme value condition of Dietrich, de Haan & Husler (2002)
#' and Husler & Li (2006)
#'
#' Tests whether the distribution is in the domain of attraction of an
#' extreme value distribution (`H0`), against the alternative that it is in
#' none, by comparing the log-spacings of the `k + 1` largest observations
#' with their limiting form. With `d(t) = log X[n-[kt],n] - log X[n-k,n]` and
#' the moment estimators `gamma_+ = M1` (Hill) and
#' `gamma_- = 1 - 1 / (2 (1 - M1^2 / M2))` of Dekkers et al. (1989), not
#' truncated at 0:
#' * `E = k integral_0^1 (d(t) / gamma_+ - (t^(-gamma_-) - 1) / gamma_- (1 -
#'   gamma_-))^2 t^eta dt` (Dietrich et al. 2002, equation 6), with
#'   `eta = 2`;
#' * `PE = k integral_0^1 (d(t) / gamma_+ + log(t))^2 t^eta dt` (equation 8,
#'   where it is called `T`; Husler & Li 2006), for `gamma >= 0`, with
#'   `eta = 0.5`, 1 or 2.
#' Logarithms are used, so the `k + 1` largest observations must be
#' positive.
#'
#' **Null distribution.** Under `H0` and the second-order conditions of
#' Dietrich et al. (Theorem 1), with `sqrt(k) A(n/k) -> 0`, both statistics
#' converge to integrals of functionals of a Brownian motion. The limit of `E`
#' depends on `gamma_-`; following Husler & Li (2006, p. 71), the quantiles
#' are taken at the moment estimate `gamma = gamma_+ + gamma_-`, interpolated
#' linearly on a grid from -0.7 to 0, with the value at 0 for `gamma >= 0` and
#' at -0.7 below it. The limit of `PE` does not depend on `gamma`. The test
#' rejects for large values. The quantiles come from a simulation of the
#' limiting variables (50000 replicates per value of `gamma`, Brownian motion
#' on 4000 points), which reproduces Tables 1 and 2 of Husler & Li (2006);
#' p-values below 5e-4 are reported as 5e-4.
#'
#' The limit of `E` in equation 7 of Dietrich et al. (2002) has a sign error
#' in its last term, which Husler & Li (2006) correct; the quantiles of
#' Dietrich et al. (Table 1) for `gamma < 0` were computed with the printed
#' sign. The tables used here follow the correction. See the reproduction
#' vignette.
#'
#' Husler & Li (2006) recommend `PE` (or the test of [drees2006()]) when the
#' estimates of `gamma` are clearly positive and `E` otherwise, with `k` at
#' most 5% of `n`.
#'
#' @param x Numeric sample; the `k + 1` largest values must be positive.
#' @param k Number of log-excesses: whole numbers from 2 to `n - 1`, or
#'   `NULL` for all of them.
#' @param statistic `"E"` or `"PE"`.
#' @param eta The exponent of the weight: 2 for `E`; 0.5, 1 or 2 for `PE`.
#'
#' @return An `evt_test` data frame with columns `k`, `statistic`, `p_value`
#'   and `gamma` (the moment estimate `gamma_+ + gamma_-`).
#'
#' @references
#' Dietrich, D., de Haan, L., & Husler, J. (2002). Testing extreme value
#'   conditions. Extremes, 5(1), 71-85. doi:10.1023/A:1020934126695
#'
#' Husler, J., & Li, D. (2006). On testing extreme value conditions.
#'   Extremes, 9(1), 69-86. doi:10.1007/s10687-006-0025-8
#'
#' @examples
#' set.seed(1)
#' x <- stats::runif(2000)^(-1)          # Pareto: gamma = 1
#' dietrich2002(x, k = c(50, 100), statistic = "PE", eta = 0.5)
#' dietrich2002(exp(stats::rexp(2000)), k = c(50, 100))
#'
#' @export
dietrich2002 <- function(x, k = NULL, statistic = c("E", "PE"),
                         eta = if (statistic == "E") 2 else 0.5) {
  statistic <- match.arg(statistic)
  allowed <- if (statistic == "E") 2 else c(0.5, 1, 2)
  if (!is.numeric(eta) || length(eta) != 1L || !eta %in% allowed) {
    stop("`eta` must be ", paste(allowed, collapse = ", "), " for the ",
         statistic, " test.", call. = FALSE)
  }
  label <- "Dietrich, de Haan & Husler (2002)"
  z <- evt_sorted(x, label, min_n = 3L)
  n <- length(z)
  k <- evt_k(k, 2L, n - 1L, label)
  if (any(z[seq_len(max(k) + 1L)] <= 0)) {
    stop("the test needs the k + 1 largest observations to be positive.",
         call. = FALSE)
  }
  lz <- log(z[seq_len(max(k) + 1L)])
  rows <- lapply(k, function(kk) {
    d <- lz[seq_len(kk)] - lz[kk + 1L]
    m1 <- mean(d)
    m2 <- mean(d^2)
    if (!(m1 > 0)) {
      stop("the k + 1 largest observations must not all be equal.",
           call. = FALSE)
    }
    gm <- 1 - 0.5 / (1 - m1^2 / m2)
    s <- ev_condition_stat(d, m1, gm, eta, statistic)
    q <- if (statistic == "E") {
      table_at(limit_tables$E, m1 + gm)
    } else {
      limit_tables$PE$q[, match(eta, limit_tables$PE$eta)]
    }
    data.frame(k = kk, statistic = s, p_value = table_pvalue(s, q),
               gamma = m1 + gm)
  })
  new_evt_test(do.call(rbind, rows),
               paste0("Dietrich, de Haan & Husler (2002), ", statistic,
                      ", eta = ", eta),
               "extreme value condition", "two.sided")
}

#' E or PE from the log-spacings d (decreasing order statistics); the
#' integrand is constant in the data on each piece [j/k, (j+1)/k)
#' @noRd
ev_condition_stat <- function(d, gp, gm, eta, statistic) {
  k <- length(d)
  a <- d / gp
  if (statistic == "PE") {
    # exact: an antiderivative of t^eta (a + log t)^2
    prim <- function(t, a) {
      if (t == 0) return(0)
      l <- a + log(t)
      e1 <- eta + 1
      t^e1 / e1 * (l^2 - 2 * l / e1 + 2 / e1^2)
    }
    s <- vapply(seq_len(k), function(j) {
      prim(j / k, a[j]) - prim((j - 1) / k, a[j])
    }, numeric(1))
    return(k * sum(s))
  }
  cfun <- function(t) if (abs(gm) < 1e-8) -log(t) else (t^(-gm) - 1) / gm
  f <- function(t, aj) (aj - cfun(t) * (1 - gm))^2 * t^eta
  # the 8-point Gauss-Legendre rule on the pieces away from 0, where the
  # integrand is smooth; adaptive quadrature on the first piece
  first <- stats::integrate(f, 0, 1 / k, aj = a[1L], rel.tol = 1e-12)$value
  if (k == 1L) return(k * first)
  gl <- gauss_legendre8
  j <- rep(1:(k - 1L), each = 8L)
  t <- (j + (rep(gl$x, times = k - 1L) + 1) / 2) / k
  w <- rep(gl$w, times = k - 1L) / (2 * k)
  k * (first + sum(w * f(t, a[j + 1L])))
}

#' Nodes and weights of the 8-point Gauss-Legendre rule on (-1, 1)
#' @noRd
gauss_legendre8 <- list(
  x = c(-0.960289856497536, -0.796666477413627, -0.525532409916329,
        -0.183434642495650, 0.183434642495650, 0.525532409916329,
        0.796666477413627, 0.960289856497536),
  w = c(0.101228536290376, 0.222381034453374, 0.313706645877887,
        0.362683783378362, 0.362683783378362, 0.313706645877887,
        0.222381034453374, 0.101228536290376)
)

#' Choulakian & Stephens (2001) goodness-of-fit tests for the generalized
#' Pareto distribution
#'
#' Tests whether a sample of excesses comes from a generalized Pareto
#' distribution with unknown shape and scale (their Case 3), with the
#' Cramer-von Mises statistic `W2` and the Anderson-Darling statistic `A2`.
#' The parameters are estimated by maximum likelihood (computed with
#' [grimshaw1993()], whose shape `k` is the paper's: `k = -gamma`), the
#' sample is transformed by the fitted distribution function, `z = F(x)`, and
#' `W2 = sum((z[i] - (2i - 1) / (2n))^2) + 1 / (12n)`,
#' `A2 = -n - (1/n) sum((2i - 1) (log z[i] + log(1 - z[n + 1 - i])))`
#' (Section 2). Large values reject.
#'
#' The p-values are obtained from the asymptotic percentage points of the
#' paper's Table 2 (to be used without modification for `n >= 25`): the
#' points are interpolated linearly in the estimated `k` (with `k = 0.5` for
#' larger estimates, as the paper says, and `k = -0.9` for smaller ones,
#' which the paper does not cover), and the p-value is interpolated linearly
#' in `log(p)` between the tabulated levels 0.5, 0.25, 0.1, 0.05, 0.025,
#' 0.01, 0.005 and 0.001. P-values beyond the table are reported at its
#' edges, 0.5 and 0.001. The test needs the maximum likelihood estimate,
#' which the paper assumes to exist.
#'
#' @param y Positive numeric vector of excesses over a threshold.
#'
#' @return An `evt_test` data frame with one row and columns `n`, `k`
#'   (Grimshaw's shape), `scale`, `W2`, `p_W2`, `A2` and `p_A2`.
#'
#' @references
#' Choulakian, V., & Stephens, M. A. (2001). Goodness-of-fit tests for the
#'   generalized Pareto distribution. Technometrics, 43(4), 478-484.
#'   doi:10.1198/00401700152672573
#'
#' @seealso [bader2018()] for the choice of a threshold by sequential
#'   testing; [wheaton_river], the example of the paper.
#'
#' @examples
#' # Section 3 of the paper: the GPD does not fit at the threshold 27.5
#' choulakian_stephens2001(wheaton_river$excess)
#'
#' @export
choulakian_stephens2001 <- function(y) {
  label <- "Choulakian & Stephens (2001)"
  fit <- grimshaw1993(y)
  if (!fit$exists) {
    stop("the maximum likelihood estimate does not exist for these data.",
         call. = FALSE)
  }
  n <- length(y)
  kk <- fit$k
  a <- fit$scale
  z <- sort(if (abs(kk) < 1e-12) 1 - exp(-y / a) else
    1 - (1 - kk * y / a)^(1 / kk))
  z <- pmin(pmax(z, 1e-15), 1 - 1e-15)
  i <- seq_len(n)
  w2 <- sum((z - (2 * i - 1) / (2 * n))^2) + 1 / (12 * n)
  a2 <- -n - mean((2 * i - 1) * (log(z) + log(1 - rev(z))))
  new_evt_test(data.frame(n = n, k = kk, scale = a,
                          W2 = w2, p_W2 = cs_pvalue(w2, kk, "W2"),
                          A2 = a2, p_A2 = cs_pvalue(a2, kk, "A2")),
               label, "generalized Pareto distribution", "two.sided")
}

#' Table 2 of Choulakian & Stephens (2001), Case 3
#' @noRd
cs_table <- list(
  k = c(-0.9, -0.5, -0.2, -0.1, 0, 0.1, 0.2, 0.3, 0.4, 0.5),
  p = c(0.5, 0.25, 0.1, 0.05, 0.025, 0.01, 0.005, 0.001),
  W2 = rbind(
    c(0.046, 0.067, 0.094, 0.115, 0.136, 0.165, 0.187, 0.239),
    c(0.049, 0.072, 0.101, 0.124, 0.147, 0.179, 0.204, 0.264),
    c(0.053, 0.078, 0.111, 0.137, 0.164, 0.200, 0.228, 0.294),
    c(0.055, 0.081, 0.116, 0.144, 0.172, 0.210, 0.240, 0.310),
    c(0.057, 0.086, 0.124, 0.153, 0.183, 0.224, 0.255, 0.330),
    c(0.059, 0.089, 0.129, 0.160, 0.192, 0.236, 0.270, 0.351),
    c(0.062, 0.094, 0.137, 0.171, 0.206, 0.254, 0.291, 0.380),
    c(0.065, 0.100, 0.147, 0.184, 0.223, 0.276, 0.317, 0.415),
    c(0.069, 0.107, 0.159, 0.201, 0.244, 0.303, 0.349, 0.458),
    c(0.074, 0.116, 0.174, 0.222, 0.271, 0.338, 0.390, 0.513)),
  A2 = rbind(
    c(0.339, 0.471, 0.641, 0.771, 0.905, 1.086, 1.226, 1.559),
    c(0.356, 0.499, 0.685, 0.830, 0.978, 1.180, 1.336, 1.707),
    c(0.376, 0.534, 0.741, 0.903, 1.069, 1.296, 1.471, 1.893),
    c(0.386, 0.550, 0.766, 0.935, 1.110, 1.348, 1.532, 1.966),
    c(0.397, 0.569, 0.796, 0.974, 1.158, 1.409, 1.603, 2.064),
    c(0.410, 0.591, 0.831, 1.020, 1.215, 1.481, 1.687, 2.176),
    c(0.426, 0.617, 0.873, 1.074, 1.283, 1.567, 1.788, 2.314),
    c(0.445, 0.649, 0.924, 1.140, 1.365, 1.672, 1.909, 2.475),
    c(0.468, 0.688, 0.985, 1.221, 1.465, 1.799, 2.058, 2.674),
    c(0.496, 0.735, 1.061, 1.321, 1.590, 1.958, 2.243, 2.922))
)

#' Percentage points at shape k (linear interpolation, clamped)
#' @noRd
cs_points <- function(k, stat) {
  tab <- cs_table[[stat]]
  kk <- min(max(k, -0.9), 0.5)
  apply(tab, 2, function(col) stats::approx(cs_table$k, col, xout = kk)$y)
}

#' p-value by interpolation in log(p) between the percentage points
#' @noRd
cs_pvalue <- function(value, k, stat) {
  q <- cs_points(k, stat)
  if (value <= q[1L]) return(0.5)
  if (value >= q[length(q)]) return(0.001)
  exp(stats::approx(q, log(cs_table$p), xout = value)$y)
}

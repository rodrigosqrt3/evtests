#' Maximum likelihood fit of the generalized Pareto distribution (Grimshaw 1993)
#'
#' Maximum likelihood estimates of the generalized Pareto distribution
#' `F(y) = 1 - (1 - k y / alpha)^(1 / k)` (`1 - exp(-y / alpha)` for
#' `k = 0`), computed by the method of Grimshaw (1993). Grimshaw's shape `k`
#' is minus the extreme value index, `gamma = -k`, and `alpha` is the scale.
#'
#' With `theta = k / alpha`, the likelihood equations reduce to a single
#' equation, `h(theta) = (1 + mean(log(1 - theta y))) mean(1 / (1 - theta y)) -
#' 1 = 0` (equation 2.2), with `k = -mean(log(1 - theta y))` (equation 2.4).
#' Every zero lies in `(theta_L, 0)` or `(0, 1 / max(y))`, with
#' `theta_L = 2 (min(y) - mean(y)) / min(y)^2` (the Theorem of p. 188). For
#' `k > 1` the likelihood is unbounded, so `k <= 1` is imposed, and the point
#' `(k, alpha) = (1, max(y))`, with log-likelihood `-n log(max(y))`, is a
#' candidate on the boundary. The estimate is the stationary point with the
#' largest likelihood if that likelihood exceeds the boundary value, and the
#' boundary point otherwise (Section 3, steps 6-8). When no stationary point
#' is a maximum, the paper concludes that there is no maximum likelihood
#' estimate.
#'
#' Grimshaw searches for at most two zeros on each side of 0 by a modified
#' Newton-Raphson method whose details and tolerances he does not give. This
#' function instead brackets every sign change of `h` on a dense grid on each
#' side (logarithmically spaced towards 0 and towards `1 / max(y)`, with the
#' paper's `eps = 1e-6 / mean(y)` as the closest distance) and refines each
#' with [stats::uniroot()], so that it does not depend on the number of
#' zeros. A zero is kept as a local maximum if the profile log-likelihood
#' (equation 2.1) is larger there than at nearby points.
#'
#' The asymptotic results quoted by Grimshaw (from Smith 1984, valid for
#' `k < 1/2`) are returned: the variance of the estimates and the
#' standard error of the upper end-point `alpha / k` when `k > 0`,
#' `alpha^2 (1 - k)^2 (1 - 2k) / (n k^4)`.
#'
#' @param y Positive numeric vector of excesses over a threshold.
#'
#' @return A list of class `evt_gpd_fit` with `gamma` (`-k`), `k`, `scale`
#'   (`alpha`), `loglik`, `theta`, `exists` (`FALSE` when there is no maximum
#'   likelihood estimate, with the estimates `NA`), `boundary` (`TRUE` when
#'   the estimate is the boundary point `k = 1`), `zeros` (a data frame of
#'   the zeros of `h` with `theta`, `k`, `alpha`, `loglik` and `maximum`),
#'   `theta_L`, `h2_0` (the limit of `h''` at 0, `mean(y^2) - 2 mean(y)^2`),
#'   `vcov` (the asymptotic covariance matrix of `(k, alpha)`, `NA` unless
#'   `k < 1/2`), `endpoint` and `endpoint_se` (`NA` unless `0 < k < 1/2`),
#'   and `n`.
#'
#' @references
#' Grimshaw, S. D. (1993). Computing maximum likelihood estimates for the
#'   generalized Pareto distribution. Technometrics, 35(2), 185-191.
#'   doi:10.1080/00401706.1993.10485040
#'
#' @seealso [drees2004()] for the fit to the excesses over a high order
#'   statistic; [nylon_fibres], the example of the paper.
#'
#' @examples
#' # Section 4 of Grimshaw (1993): k = .117698, alpha = .283040
#' fit <- grimshaw1993(nylon_fibres$excess)
#' fit
#' fit$zeros
#'
#' @export
grimshaw1993 <- function(y) {
  if (!is.numeric(y) || anyNA(y) || any(!is.finite(y))) {
    stop("`y` must be a numeric vector without missing or infinite values.",
         call. = FALSE)
  }
  if (any(y < 0) || !any(y > 0)) {
    stop("`y` must be non-negative excesses, not all zero.", call. = FALSE)
  }
  n <- length(y)
  if (n < 2L) stop("Grimshaw (1993) needs at least 2 excesses.", call. = FALSE)
  ybar <- mean(y)
  ymax <- max(y)
  ymin <- min(y)
  eps <- 1e-6 / ybar
  theta_l <- 2 * (ymin - ybar) / ymin^2
  theta_u <- 1 / ymax - eps
  h2_0 <- mean(y^2) - 2 * ybar^2

  ty <- function(th) {
    m <- outer(th, y)
    m[m >= 1] <- NA_real_          # outside the support: no log of a negative
    m
  }
  h <- function(th) {                                   # equation 2.2
    m <- ty(th)
    (1 + rowMeans(log1p(-m))) * rowMeans(1 / (1 - m)) - 1
  }
  prof <- function(th) {                                # equation 2.1
    s <- rowSums(log1p(-ty(th)))
    -n - s - n * log(-s / (n * th))
  }

  grid_side <- function(a, b, near_a) {
    # points from a to b, dense (log spacing) near the end given
    u <- exp(seq(log(eps), log(abs(b - a)), length.out = 1500L))
    pts <- if (near_a) a + sign(b - a) * u else b - sign(b - a) * u
    sort(unique(c(a, pts, b)))
  }
  zeros <- numeric(0)
  if (is.finite(theta_l) && theta_l < -eps) {
    neg <- sort(unique(c(grid_side(theta_l, -eps, TRUE),
                         grid_side(theta_l, -eps, FALSE))))
    zeros <- c(zeros, gpd_roots(h, neg))
  }
  if (theta_u > eps) {
    pos <- sort(unique(c(grid_side(eps, theta_u, TRUE),
                         grid_side(eps, theta_u, FALSE))))
    zeros <- c(zeros, gpd_roots(h, pos))
  }

  zk <- vapply(zeros, function(t) -mean(log1p(-t * y)), numeric(1))
  za <- zk / zeros
  zl <- prof(zeros)
  is_max <- vapply(seq_along(zeros), function(i) {
    d <- max(abs(zeros[i]) * 1e-4, eps)
    lo <- zeros[i] - d
    hi <- min(zeros[i] + d, 1 / ymax - eps / 10)
    zl[i] >= prof(lo) && zl[i] >= prof(hi)
  }, logical(1))
  zero_tab <- data.frame(theta = zeros, k = zk, alpha = za, loglik = zl,
                         maximum = is_max)

  cand <- which(is_max & zk <= 1)
  l_bound <- -n * log(ymax)
  out <- list(gamma = NA_real_, k = NA_real_, scale = NA_real_,
              loglik = NA_real_, theta = NA_real_, exists = FALSE,
              boundary = FALSE, zeros = zero_tab, theta_L = theta_l,
              h2_0 = h2_0, vcov = matrix(NA_real_, 2, 2), endpoint = NA_real_,
              endpoint_se = NA_real_, n = n)
  if (length(cand)) {
    best <- cand[which.max(zl[cand])]
    out$exists <- TRUE
    if (zl[best] > l_bound) {
      out[c("k", "scale", "loglik", "theta")] <-
        list(zk[best], za[best], zl[best], zeros[best])
    } else {
      out[c("k", "scale", "loglik", "theta", "boundary")] <-
        list(1, ymax, l_bound, 1 / ymax, TRUE)
    }
    out$gamma <- -out$k
    kk <- out$k
    if (kk < 0.5) {
      a <- out$scale
      out$vcov <- matrix(c((1 - kk)^2, a * (1 - kk), a * (1 - kk),
                           2 * a^2 * (1 - kk)), 2, 2) / n
      if (kk > 0) {
        out$endpoint <- a / kk
        out$endpoint_se <- sqrt(a^2 * (1 - kk)^2 * (1 - 2 * kk) / (n * kk^4))
      }
    }
  }
  dimnames(out$vcov) <- list(c("k", "alpha"), c("k", "alpha"))
  class(out) <- "evt_gpd_fit"
  out
}

#' Zeros of a function bracketed on a grid
#' @noRd
gpd_roots <- function(f, grid) {
  v <- f(grid)
  ok <- is.finite(v)
  grid <- grid[ok]
  v <- v[ok]
  i <- which(v[-1L] * v[-length(v)] < 0)
  vapply(i, function(j) {
    stats::uniroot(f, c(grid[j], grid[j + 1L]), tol = 1e-12)$root
  }, numeric(1))
}

#' Print a generalized Pareto fit
#'
#' @param x An `evt_gpd_fit` object.
#' @param digits Number of significant digits.
#' @param ... Unused.
#' @return `x`, invisibly.
#' @export
print.evt_gpd_fit <- function(x, digits = 6L, ...) {
  cat("<Generalized Pareto maximum likelihood fit, Grimshaw (1993)>\n")
  cat("  ", x$n, " excesses\n", sep = "")
  if (!x$exists) {
    cat("  no maximum likelihood estimate exists\n")
    return(invisible(x))
  }
  cat("  extreme value index gamma = ", format(x$gamma, digits = digits),
      " (Grimshaw's k = ", format(x$k, digits = digits), ")\n", sep = "")
  cat("  scale = ", format(x$scale, digits = digits), "\n", sep = "")
  cat("  log-likelihood = ", format(x$loglik, digits = digits),
      if (x$boundary) " (boundary point k = 1)", "\n", sep = "")
  invisible(x)
}

#' Drees & Kaufmann (1998) choice of the number of excesses for the Hill
#' estimator
#'
#' Estimates the number of upper order statistics that minimises the
#' asymptotic mean squared error of the Hill estimator, for a heavy tail
#' (`gamma > 0`), by the sequential procedure of Drees & Kaufmann (1998).
#' With `b(i)` the Hill estimate from the `i` largest observations, the
#' smallest of which is the threshold (so `i - 1` log-excesses), the stopping
#' time is `kbar(r) = min{k : max_{2 <= i <= k} sqrt(i) |b(i) - b(k)| > r}`
#' (equation 4), and
#' `k_opt = [(2 rho + 1)^(-1/rho) (2 g^2 rho)^(1 / (2 rho + 1))
#' (kbar(r^xi) / kbar(r)^xi)^(1 / (1 - xi))]` (Theorem 1), where `g` is an
#' initial estimate of `gamma`.
#'
#' The recommended choices of Section 3 are used: `2.5 g n^0.25` for `r`,
#' `xi = 0.7`, `g` the Hill estimate from the `2 sqrt(n+)` largest
#' observations (`n+` the number of positive observations; rounded down),
#' and `rho = 1` unless `rho = "estimate"`. When `kbar(r)` is not defined, `r`
#' is replaced by `0.9 r` until it is. With `rho = "estimate"`, `rho` is
#' estimated by `rho^(1)` of Theorem 1 with `lambda = 0.6` at the threshold
#' `r^xi`, and set to 1 when that estimate is not defined or not positive;
#' the paper calls the estimator used in its Tables 5 and 6 "modified"
#' without defining the modification, so those tables are not a target.
#' Only the positive observations are used.
#'
#' **Convention for `k`.** The paper's `k` counts the threshold; the returned
#' `k` is the number of excesses (the paper's `k` minus one), so that
#' `hill1975(x, k)` gives the adaptive estimate.
#'
#' @param x Numeric sample.
#' @param rho `1` (default, the paper's recommendation), another positive
#'   number, or `"estimate"`.
#' @param xi,r_factor The constants `xi` and 2.5 of `2.5 g n^0.25` for `r`.
#'
#' @return A list of class `evt_k_choice` with `k` (number of excesses, `NA`
#'   when the paper's `k_opt` is below 2), `gamma` (the Hill estimate at
#'   `k`), `rho`, `kbar` (the stopping times at `r` and `r^xi`, in the
#'   paper's convention), `r` and `gamma_initial`.
#'
#' @references
#' Drees, H., & Kaufmann, E. (1998). Selecting the optimal sample fraction
#'   in univariate extreme value estimation. Stochastic Processes and their
#'   Applications, 75(2), 149-172. doi:10.1016/S0304-4149(98)00017-9
#'
#' @examples
#' set.seed(1)
#' x <- (-log(stats::runif(1000)))^(-0.5)    # Frechet, gamma = 0.5
#' drees_kaufmann1998(x)
#'
#' @export
drees_kaufmann1998 <- function(x, rho = 1, xi = 0.7, r_factor = 2.5) {
  if (!(identical(rho, "estimate") ||
        (is.numeric(rho) && length(rho) == 1L && rho > 0))) {
    stop("`rho` must be a positive number or \"estimate\".", call. = FALSE)
  }
  if (!is.numeric(xi) || length(xi) != 1L || !(xi > 0 && xi < 1)) {
    stop("`xi` must be in (0, 1).", call. = FALSE)
  }
  label <- "Drees & Kaufmann (1998)"
  z <- evt_sorted(x, label, min_n = 3L)
  n <- length(z)
  zp <- z[z > 0]
  np <- length(zp)
  if (np < 4L) {
    stop(label, " needs at least 4 positive observations.", call. = FALSE)
  }
  lz <- log(zp)
  # b[i] = Hill from the i largest (threshold zp[i]), i = 2, ..., np
  cs <- cumsum(lz)
  b <- c(NA_real_, cs[-np] / seq_len(np - 1L) - lz[-1L])
  m0 <- max(2L, floor(2 * sqrt(np)))
  g0 <- b[min(m0, np)]
  if (!(g0 > 0)) stop("the initial Hill estimate is not positive.",
                      call. = FALSE)
  kbar <- function(r) {
    repeat {
      for (k in 2:np) {
        i <- 2:k
        if (max(sqrt(i) * abs(b[i] - b[k])) > r) return(c(k = k, r = r))
      }
      r <- 0.9 * r
    }
  }
  r <- r_factor * g0 * n^0.25
  k1 <- kbar(r)
  k2 <- kbar(r^xi)
  rho_used <- if (identical(rho, "estimate")) {
    dk_rho1(b, k2[["k"]], 0.6)
  } else rho
  kopt <- floor((2 * rho_used + 1)^(-1 / rho_used) *
                  (2 * g0^2 * rho_used)^(1 / (2 * rho_used + 1)) *
                  (k2[["k"]] / k1[["k"]]^xi)^(1 / (1 - xi)))
  kopt <- min(kopt, np)
  k <- if (kopt >= 2L) kopt - 1L else NA_integer_
  structure(list(k = k,
                 gamma = if (is.na(k)) NA_real_ else b[kopt],
                 rho = rho_used,
                 kbar = c(r = k1[["k"]], r_xi = k2[["k"]]),
                 r = c(r = k1[["r"]], r_xi = k2[["r"]]),
                 gamma_initial = g0, method = label),
            class = "evt_k_choice")
}

#' rho^(1) of Theorem 1 of Drees & Kaufmann (1998); 1 when not defined
#' @noRd
dk_rho1 <- function(b, kbar, lambda) {
  kl <- floor(lambda * kbar)
  if (kl < 3L) return(1)
  num <- max(sqrt(2:kl) * abs(b[2:kl] - b[kl]))
  den <- max(sqrt(2:kbar) * abs(b[2:kbar] - b[kbar]))
  est <- log(num / den) / log(lambda) - 0.5
  if (!is.finite(est) || est <= 0) 1 else est
}

#' Print a choice of k
#'
#' @param x An `evt_k_choice` object.
#' @param ... Unused.
#' @return `x`, invisibly.
#' @export
print.evt_k_choice <- function(x, ...) {
  cat("<", x$method, ">\n", sep = "")
  cat("  selected k (excesses): ", x$k, "\n", sep = "")
  if (!is.null(x$gamma)) {
    cat("  gamma at k: ", format(x$gamma, digits = 4), "\n", sep = "")
  }
  invisible(x)
}

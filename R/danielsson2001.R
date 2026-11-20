#' Danielsson, de Haan, Peng & de Vries (2001) bootstrap choice of the number
#' of excesses for the Hill estimator
#'
#' Chooses the number `k` of excesses that minimises the asymptotic mean
#' squared error of the Hill estimator of a heavy tail (`gamma > 0`) by the
#' double subsample bootstrap of Danielsson et al. (2001). With
#' `gamma(k)` the Hill estimator and `M(k)` the mean squared log-excess over
#' `X[n-k,n]`, the bootstrap criterion for subsamples of size `n1` drawn with
#' replacement is `Q(n1, k1) = E((M*(k1) - 2 gamma*(k1)^2)^2)` (p. 230),
#' estimated by the average over `B` resamples; `k1` minimises it. The same is
#' done for `n2 = n1^2 / n`, giving `k2`, and (Corollary 7)
#' `k0 = (k1^2 / k2) ((log k1)^2 / (2 log n1 - log k1)^2)^((log n1 - log k1) /
#' log n1)`, with `rho = log(k1) / (2 log(k1) - 2 log(n1))` (equation 3.9).
#' `n1` is chosen on a grid to minimise `R(n1) = Q(n1, k1)^2 / Q(n2, k2)`
#' (pp. 232-233).
#'
#' The paper uses, for `n = 2000`, the grid `n1 = 600, 700, ..., 1700` with
#' `B = 1000` resamples, and for `n = 20000` the grid `2000, 3000, ..., 15000`
#' with `B = 500`; the default grid here is the first one scaled to `n`
#' (`n1 = 0.3 n, 0.35 n, ..., 0.85 n`, rounded). `n2` is rounded to the
#' nearest integer, and `k0` too (the paper prints no rounding). The theory
#' requires `rho < 0` and `n1 = O(n^(1 - e))` with `0 < e < 1/2`. Only the
#' positive observations of each resample enter the Hill estimator. The
#' resamples are drawn from R's random number generator, so set a seed for
#' reproducible results.
#'
#' @param x Numeric sample.
#' @param n1 Grid of subsample sizes (whole numbers between 2 and `n - 1`);
#'   by default `round(n * seq(0.3, 0.85, by = 0.05))`.
#' @param B Number of bootstrap resamples for each subsample size.
#'
#' @return A list of class `evt_k_choice` with `k` (`k0`, the number of
#'   excesses), `gamma` (the Hill estimate at `k`), `rho`, `n1` (the selected
#'   subsample size), `k1`, `k2`, and `grid`, a data frame with `n1`, `n2`,
#'   `k1`, `k2` and `R` for every subsample size.
#'
#' @references
#' Danielsson, J., de Haan, L., Peng, L., & de Vries, C. G. (2001). Using a
#'   bootstrap method to choose the sample fraction in tail index
#'   estimation. Journal of Multivariate Analysis, 76(2), 226-248.
#'   doi:10.1006/jmva.2000.1903
#'
#' @examples
#' set.seed(1)
#' x <- (-log(stats::runif(2000)))^(-0.25)   # Frechet, gamma = 0.25
#' danielsson2001(x, B = 200)
#'
#' @export
danielsson2001 <- function(x, n1 = NULL, B = 1000) {
  label <- "Danielsson et al. (2001)"
  z <- evt_sorted(x, label, min_n = 10L)
  n <- length(z)
  if (is.null(n1)) n1 <- unique(round(n * seq(0.3, 0.85, by = 0.05)))
  if (!is.numeric(n1) || anyNA(n1) || any(n1 != round(n1)) || any(n1 < 2) ||
      any(n1 >= n)) {
    stop("`n1` must be whole numbers between 2 and n - 1.", call. = FALSE)
  }
  if (!is.numeric(B) || length(B) != 1L || B < 1 || B != round(B)) {
    stop("`B` must be a positive whole number.", call. = FALSE)
  }
  grid <- do.call(rbind, lapply(n1, function(m1) {
    m2 <- max(2, round(m1^2 / n))
    a <- dan_q(z, m1, B)
    b <- dan_q(z, m2, B)
    data.frame(n1 = m1, n2 = m2, k1 = a$k, k2 = b$k,
               R = a$q^2 / b$q)
  }))
  ok <- which(is.finite(grid$R))
  if (!length(ok)) stop("the bootstrap criterion could not be evaluated.",
                        call. = FALSE)
  best <- grid[ok[which.min(grid$R[ok])], ]
  k1 <- best$k1
  k2 <- best$k2
  ln1 <- log(best$n1)
  k0 <- (k1^2 / k2) *
    ((log(k1))^2 / (2 * ln1 - log(k1))^2)^((ln1 - log(k1)) / ln1)
  k0 <- min(max(1L, round(k0)), n - 1L)
  zp <- z[z > 0]
  g <- if (k0 < length(zp)) mean(log(zp[seq_len(k0)])) - log(zp[k0 + 1L]) else
    NA_real_
  structure(list(k = as.integer(k0), gamma = g,
                 rho = log(k1) / (2 * log(k1) - 2 * ln1),
                 n1 = best$n1, k1 = k1, k2 = k2, grid = grid,
                 method = label),
            class = "evt_k_choice")
}

#' Bootstrap criterion Q(m, k) for k = 1, ..., m - 1 and its minimiser
#' @noRd
dan_q <- function(z, m, B) {
  acc <- numeric(m - 1L)
  cnt <- numeric(m - 1L)
  for (b in seq_len(B)) {
    s <- sort(sample(z, m, replace = TRUE), decreasing = TRUE)
    s <- s[s > 0]
    kk <- length(s) - 1L
    if (kk < 1L) next
    l <- log(s)
    cs <- cumsum(l)[seq_len(kk)]
    cs2 <- cumsum(l^2)[seq_len(kk)]
    k <- seq_len(kk)
    thr <- l[k + 1L]
    g <- cs / k - thr
    mm <- (cs2 - 2 * thr * cs) / k + thr^2
    acc[k] <- acc[k] + (mm - 2 * g^2)^2
    cnt[k] <- cnt[k] + 1
  }
  q <- ifelse(cnt == B, acc / cnt, NA_real_)
  if (all(is.na(q))) return(list(k = NA_integer_, q = NA_real_))
  k <- which.min(q)
  list(k = k, q = q[k])
}

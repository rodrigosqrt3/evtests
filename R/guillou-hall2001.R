#' Guillou & Hall (2001) diagnostic for the number of excesses of the Hill
#' estimator
#'
#' Chooses the number `k` of upper order statistics for the Hill estimator of
#' a heavy tail (`gamma > 0`) with the diagnostic of Guillou & Hall (2001).
#' With the scaled log-spacings `Y[i] = i (log X[n-i+1,n] - log X[n-i,n])`,
#' `i = 1, ..., k`, the Hill estimate `gamma_k = mean(Y[1:k])` (`k` excesses
#' over `X[n-k,n]`) and the antisymmetric weights
#' `w[i] = sign(k - 2i + 1) |k - 2i + 1|^p`,
#' `T_k = sqrt((2p + 1) / k^(2p + 1)) sum(w[i] Y[i]) / gamma_k` (equation
#' 2.3), which is approximately standard normal while the bias of the Hill
#' estimator is negligible, and its smoothed version
#' `Q_k = sqrt(sum_{j=-l}^{l} T_{k+j}^2 / (2l + 1))`, `l = floor(k / 2)`
#' (equation 2.6). The selected `k` is where the bias becomes significant:
#' with `rule = "first"`, the least `k >= kmin` with `Q_k >= crit` (Section
#' 4.2); with `rule = "all"`, the least `k` such that `Q_l > crit` for every
#' `l >= k` examined (Section 2.2, stated there for `T`). The paper recommends
#' `Q`, `p = 1` and `crit = 1.25` (or 1.5), and the paper is written for the
#' lower tail of the reciprocals, which is the upper tail used here.
#'
#' The two sections of the paper state different rules and do not give the
#' smallest `k` at which the search starts. With the first crossing from
#' `k = 2`, the mean selected `k` of Table 1 is roughly reproduced for the
#' Frechet and Burr distributions, but the relative errors `R(c)` of the
#' table (close to 1) are not: about 20% of the samples cross at small `k`.
#' See the reproduction vignette.
#'
#' @param x Numeric sample; the observations used must be positive.
#' @param crit The critical value.
#' @param rule `"first"` or `"all"`, see above.
#' @param kmin Smallest `k` for the rule `"first"`.
#' @param p The exponent of the weights.
#' @param statistic `"Q"` (default) or `"T"` (then `|T_k|` is compared with
#'   `crit`).
#'
#' @return A list of class `evt_k_choice` with `k` (`NA` when the critical
#'   value is never reached), `gamma` (the Hill estimate at `k`), and `path`,
#'   a data frame with columns `k`, `gamma`, `T` and `Q` (`Q` is `NA` where
#'   the window passes the end of the sample).
#'
#' @references
#' Guillou, A., & Hall, P. (2001). A diagnostic for selecting the threshold
#'   in extreme value analysis. Journal of the Royal Statistical Society,
#'   Series B, 63(2), 293-305. doi:10.1111/1467-9868.00286
#'
#' @examples
#' set.seed(1)
#' x <- (-log(stats::runif(500)))^(-1)       # Frechet, gamma = 1
#' guillou_hall2001(x)
#'
#' @export
guillou_hall2001 <- function(x, crit = 1.25, rule = c("first", "all"),
                             kmin = 2L, p = 1, statistic = c("Q", "T")) {
  rule <- match.arg(rule)
  statistic <- match.arg(statistic)
  label <- "Guillou & Hall (2001)"
  if (!is.numeric(crit) || length(crit) != 1L || !(crit > 0)) {
    stop("`crit` must be a positive number.", call. = FALSE)
  }
  z <- evt_sorted(x, label, min_n = 4L)
  zp <- z[z > 0]
  m <- length(zp) - 1L                       # largest k with positive data
  if (m < 3L) stop(label, " needs at least 4 positive observations.",
                   call. = FALSE)
  lz <- log(zp)
  y <- seq_len(m) * (lz[seq_len(m)] - lz[seq_len(m) + 1L])
  ks <- seq_len(m)
  gam <- cumsum(y) / ks
  tk <- vapply(ks, function(k) {
    i <- seq_len(k)
    w <- sign(k - 2 * i + 1) * abs(k - 2 * i + 1)^p
    sqrt((2 * p + 1) / k^(2 * p + 1)) * sum(w * y[i]) / gam[k]
  }, numeric(1))
  qk <- vapply(ks, function(k) {
    l <- k %/% 2L
    if (k + l > m || k - l < 1L) return(NA_real_)
    sqrt(mean(tk[(k - l):(k + l)]^2))
  }, numeric(1))
  path <- data.frame(k = ks, gamma = gam, T = tk, Q = qk)
  s <- if (statistic == "Q") qk else abs(tk)
  sel <- NA_integer_
  if (rule == "first") {
    hit <- which(ks >= kmin & !is.na(s) & s >= crit)
    if (length(hit)) sel <- ks[hit[1L]]
  } else {
    ok <- which(!is.na(s))
    below <- ok[s[ok] <= crit]
    if (length(below) && max(below) < max(ok)) sel <- ks[max(below) + 1L]
  }
  structure(list(k = sel, gamma = if (is.na(sel)) NA_real_ else gam[sel],
                 path = path, method = label),
            class = "evt_k_choice")
}

#' Bader, Yan & Zhang (2018) threshold choice by ordered goodness-of-fit
#' tests with ForwardStop
#'
#' Chooses the threshold above which the excesses follow a generalized Pareto
#' distribution by testing a sequence of candidate thresholds, from the lowest
#' to the highest, and controlling the false discovery rate with the
#' ForwardStop rule (equation 2): with `p[1], ..., p[l]` the goodness-of-fit
#' p-values of the thresholds in increasing order,
#' `kF = max{j : -(1/j) sum_{i <= j} log(1 - p[i]) <= alpha}`; the
#' hypotheses `1, ..., kF` are rejected and the selected threshold is the
#' next one, `u[kF + 1]`, or the lowest threshold when nothing is rejected
#' and the highest when everything is (p. 320).
#'
#' The candidate thresholds are the order statistics `X[n-k,n]` for the
#' numbers of excesses `k` given (lowest threshold = largest `k`). The
#' goodness-of-fit tests are those of [choulakian_stephens2001()]
#' (Anderson-Darling by default, as recommended by Bader et al.). Bader et al.
#' compute the p-values from a table of their own, simulated for
#' `n = 1000` and available only in their software; this function uses the
#' asymptotic table of Choulakian & Stephens (2001), on which the paper's
#' description of the null distributions rests. The p-values of successive
#' thresholds are dependent, which the paper studies by simulation.
#'
#' @param x Numeric sample.
#' @param k Candidate numbers of excesses (at least 10 each), or `NULL` for
#'   the 20 values `round(seq(n / 2, max(25, n / 50), length.out = 20))`.
#' @param alpha Level of the ForwardStop rule.
#' @param statistic `"A2"` or `"W2"`.
#'
#' @return A list of class `evt_k_choice` with `k` (the number of excesses at
#'   the selected threshold), `threshold`, `kF` (the number of rejected
#'   thresholds) and `path`, a data frame with `k`, `threshold`, `p_value`
#'   and the ForwardStop statistic `forward_stop`, ordered from the lowest
#'   threshold.
#'
#' @references
#' Bader, B., Yan, J., & Zhang, X. (2018). Automated threshold selection for
#'   extreme value analysis via ordered goodness-of-fit tests with
#'   adjustment for false discovery rate. The Annals of Applied Statistics,
#'   12(1), 310-329. doi:10.1214/17-AOAS1092
#'
#' @examples
#' set.seed(1)
#' x <- c(stats::rexp(500, 2), 2 + stats::rexp(200))   # GPD only in the tail
#' bader2018(x)
#'
#' @export
bader2018 <- function(x, k = NULL, alpha = 0.05,
                      statistic = c("A2", "W2")) {
  statistic <- match.arg(statistic)
  check_alpha(alpha)
  label <- "Bader, Yan & Zhang (2018)"
  z <- evt_sorted(x, label, min_n = 11L)
  n <- length(z)
  if (is.null(k)) {
    k <- unique(round(seq(n / 2, max(25, n / 50), length.out = 20)))
  }
  k <- sort(evt_k(k, 10L, n - 1L, label), decreasing = TRUE)
  p <- vapply(k, function(kk) {
    y <- evt_excess(z, kk)
    r <- tryCatch(choulakian_stephens2001(y), error = function(e) NULL)
    if (is.null(r)) NA_real_ else r[[paste0("p_", statistic)]]
  }, numeric(1))
  if (anyNA(p)) {
    stop("the maximum likelihood fit does not exist for k = ",
         paste(k[is.na(p)], collapse = ", "), "; remove these candidates.",
         call. = FALSE)
  }
  fs <- forward_stop(p)
  kf <- max(c(0L, which(fs <= alpha)))
  sel <- if (kf == length(k)) length(k) else kf + 1L
  path <- data.frame(k = k, threshold = z[k + 1L], p_value = p,
                     forward_stop = fs)
  structure(list(k = k[sel], threshold = z[k[sel] + 1L], kF = kf,
                 path = path, method = label),
            class = "evt_k_choice")
}

#' ForwardStop statistics -(1/j) sum_{i <= j} log(1 - p_i), j = 1, ..., l
#' @noRd
forward_stop <- function(p) -cumsum(log1p(-p)) / seq_along(p)

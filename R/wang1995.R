#' Wang (1995) choice of the number of upper order statistics
#'
#' Chooses how many of the largest observations follow the Gumbel-domain
#' limit, by applying the Hasofer & Wang (1992) statistic `W` to the
#' `k + 1` largest observations for `k = 2, 3, ...` and stopping at the first
#' rejection: the selected `k` is the last one before the first rejection, so
#' that the fit is acceptable for every smaller `k` (Section 2, p. 1056).
#' Wang counts the threshold among the top order statistics; his selected
#' number of observations is `k + 1` here.
#'
#' `W` is tested two-sided, with the limiting null of [hasofer_wang1992()]
#' (Section 3.2; Wang multiplies `W` by the number of observations, which does
#' not change the test). The paper does not say how `alpha` is shared between
#' the two tails; the published examples are reproduced with equal tails
#' (`alpha / 2` each), which is what this function does. For the Frechet
#' domain the paper applies the procedure to the logarithms of the
#' observations (`transform = "log"`), and for the Weibull domain with a known
#' end-point `q` to `-log(q - x)` (`transform = "endpoint"`). Wang recommends
#' `alpha` of at least 10% for samples of up to about 500 observations. The
#' variant with an estimated end-point and the R-squared criterion of
#' Section 3.1 are not implemented: the first needs an end-point estimator that
#' the paper takes from other work, and the critical values of the second are
#' not given.
#'
#' @param x Numeric sample.
#' @param alpha Level of each test, in (0, 1).
#' @param transform `"none"`, `"log"` (Frechet domain) or `"endpoint"`
#'   (Weibull domain with known end-point `q`).
#' @param q The end-point, for `transform = "endpoint"`; it must exceed the
#'   sample maximum.
#' @param kmax Largest `k` examined; by default `n - 1`.
#' @param nsim Number of simulations of the null distribution of `W`.
#'
#' @return An `evt_test` data frame with the path (`k`, `statistic` and
#'   two-sided `p_value` for `k = 2, ..., kmax`) and the attribute
#'   `selected`, the chosen `k` (`kmax` when no test rejects).
#'
#' @references
#' Wang, J. Z. (1995). Selection of the k largest order statistics for the
#'   domain of attraction of the Gumbel distribution. Journal of the
#'   American Statistical Association, 90(431), 1055-1061.
#'   doi:10.1080/01621459.1995.10476607
#'
#' @seealso [bushfires] and [mice_survival], the examples of the paper.
#'
#' @examples
#' # Example 1 of the paper, Frechet domain: 27 top observations, so k = 26
#' w <- wang1995(bushfires$area, alpha = 0.1, transform = "log", nsim = 2e4)
#' attr(w, "selected")
#'
#' @export
wang1995 <- function(x, alpha = 0.1,
                     transform = c("none", "log", "endpoint"), q = NULL,
                     kmax = NULL, nsim = 1e5) {
  check_alpha(alpha)
  transform <- match.arg(transform)
  y <- switch(transform,
    none = x,
    log = {
      if (any(x <= 0)) stop("`transform = \"log\"` needs positive data.",
                            call. = FALSE)
      log(x)
    },
    endpoint = {
      if (is.null(q) || !is.numeric(q) || length(q) != 1L || !(q > max(x))) {
        stop("`q` must be a single number larger than the sample maximum.",
             call. = FALSE)
      }
      -log(q - x)
    })
  z <- evt_sorted(y, "Wang (1995)", min_n = 3L)
  n <- length(z)
  if (is.null(kmax)) kmax <- n - 1L
  k <- evt_k(seq.int(2L, kmax), 2L, n - 1L, "Wang (1995)")
  path <- hasofer_wang1992(z, k = k, nsim = nsim)
  rej <- which(path$p_value < alpha)
  selected <- if (length(rej)) path$k[rej[1L]] - 1L else max(k)
  out <- new_evt_test(as.data.frame(path), "Wang (1995)",
                      "gamma = 0 (Gumbel domain) for each k", "two.sided")
  attr(out, "selected") <- selected
  out
}

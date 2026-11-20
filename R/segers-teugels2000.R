#' Segers & Teugels (2000) test of the Gumbel domain with Galton's ratio
#'
#' Tests `H0: gamma = 0` (the Gumbel domain) with Galton's ratio of the three
#' largest observations of each of `B` subsamples. In subsample `i`,
#' `xi[i] = (X(1) - X(3)) / (X(2) - X(3))` with `X(1) >= X(2) >= X(3)` its
#' three largest values; under `H0` (in the limit) `xi[i]` has the
#' distribution `1 - 2 / (1 + x)`, `x > 1`, free of any parameter (Corollary
#' 2.2). With the score `T(x) = 1 - 6 x / (1 + x)^2` (equation 3.5), which has
#' mean 0 and variance 1/5 under `H0`:
#' * two-sided: `S = (5 / B) (sum T(xi[i]))^2`, compared with chi-squared on
#'   one degree of freedom;
#' * one-sided: `Z = sqrt(5 / B) sum T(xi[i])`, compared with the standard
#'   normal; large `Z` points to `gamma > 0` and small `Z` to `gamma < 0`.
#' Both are asymptotic as the number of subsamples grows.
#'
#' The sample is divided, in the order given, into consecutive subsamples of
#' `block_size` observations (the paper's simulations use 50 subsamples of
#' 20); observations left over at the end are not used. The paper does not
#' say how to form the subsamples or how to handle ties. The statistics are
#' location and scale invariant.
#'
#' @param x Numeric sample, in its natural order (for instance time).
#' @param block_size Number of observations in each subsample, at least 3.
#' @param alternative `"two.sided"` (the statistic `S`), `"positive"` or
#'   `"negative"` (the statistic `Z`).
#'
#' @return An `evt_test` data frame with one row and columns `blocks`,
#'   `statistic` (`S` or `Z`) and `p_value`.
#'
#' @references
#' Segers, J., & Teugels, J. (2000). Testing the Gumbel hypothesis by
#'   Galton's ratio. Extremes, 3(3), 291-303. doi:10.1023/A:1011411212045
#'
#' @examples
#' set.seed(1)
#' segers_teugels2000(stats::rexp(1000))
#' segers_teugels2000(stats::runif(1000)^(-0.5), alternative = "positive")
#'
#' @export
segers_teugels2000 <- function(x, block_size = 20L,
                               alternative = c("two.sided", "positive",
                                               "negative")) {
  alternative <- evt_alternative(alternative)
  label <- "Segers & Teugels (2000)"
  if (!is.numeric(x) || anyNA(x) || any(!is.finite(x))) {
    stop("`x` must be a numeric vector without missing or infinite values.",
         call. = FALSE)
  }
  if (!is.numeric(block_size) || length(block_size) != 1L ||
      block_size < 3 || block_size != round(block_size)) {
    stop("`block_size` must be a whole number of at least 3.", call. = FALSE)
  }
  nb <- length(x) %/% block_size
  if (nb < 2L) stop(label, " needs at least two subsamples.", call. = FALSE)
  xi <- vapply(seq_len(nb), function(i) {
    b <- sort(x[(i - 1) * block_size + seq_len(block_size)],
              decreasing = TRUE)[1:3]
    if (!(b[2L] > b[3L])) {
      stop("the second and third largest values of a subsample are tied.",
           call. = FALSE)
    }
    (b[1L] - b[3L]) / (b[2L] - b[3L])
  }, numeric(1))
  tt <- st_score(xi)
  if (alternative == "two.sided") {
    s <- (5 / nb) * sum(tt)^2
    p <- 1 - stats::pchisq(s, 1)
  } else {
    s <- sqrt(5 / nb) * sum(tt)
    p <- if (alternative == "positive") 1 - stats::pnorm(s) else
      stats::pnorm(s)
  }
  new_evt_test(data.frame(blocks = nb, statistic = s, p_value = p), label,
               "gamma = 0 (Gumbel domain)", alternative)
}

#' Score of equation 3.5
#' @noRd
st_score <- function(x) 1 - 6 * x / (1 + x)^2

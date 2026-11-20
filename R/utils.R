#' Validate a sample and return it sorted in decreasing order
#'
#' @param x Numeric sample.
#' @param label Name of the method, for error messages.
#' @param min_n Minimum sample size.
#' @param positive If `TRUE`, the observations used must be positive (checked
#'   by the caller on the order statistics actually used).
#' @return The finite values of `x`, sorted in decreasing order.
#' @noRd
evt_sorted <- function(x, label, min_n = 2L) {
  if (!is.numeric(x)) stop("`x` must be a numeric vector.", call. = FALSE)
  if (anyNA(x) || any(!is.finite(x))) {
    stop("`x` must not contain missing or infinite values.", call. = FALSE)
  }
  n <- length(x)
  if (n < min_n) {
    stop(label, " needs at least ", min_n, " observations; `x` has ", n, ".",
         call. = FALSE)
  }
  sort(x, decreasing = TRUE)
}

#' Validate the number of upper order statistics
#'
#' @param k `NULL` or a vector of whole numbers.
#' @param kmin,kmax Admissible range.
#' @param label Name of the method, for error messages.
#' @return Sorted unique integer vector; `kmin:kmax` when `k` is `NULL`.
#' @noRd
evt_k <- function(k, kmin, kmax, label) {
  if (kmax < kmin) {
    stop(label, " needs a larger sample: the admissible k range is empty.",
         call. = FALSE)
  }
  if (is.null(k)) return(seq.int(kmin, kmax))
  if (!is.numeric(k) || !length(k) || anyNA(k) || any(k != round(k))) {
    stop("`k` must be whole numbers.", call. = FALSE)
  }
  if (any(k < kmin) || any(k > kmax)) {
    stop("`k` must be between ", kmin, " and ", kmax, " for ", label,
         " with this sample.", call. = FALSE)
  }
  sort(unique(as.integer(k)))
}

#' Validate a significance level
#' @noRd
check_alpha <- function(alpha) {
  if (!is.numeric(alpha) || length(alpha) != 1L || is.na(alpha) ||
      alpha <= 0 || alpha >= 1) {
    stop("`alpha` must be a single number in (0, 1).", call. = FALSE)
  }
  invisible(alpha)
}

#' Construct an estimate path over k
#'
#' @param df Data frame with a column `k` and one row per value of k.
#' @param method Name of the method.
#' @param param Name of the main column (plotted by default).
#' @noRd
new_evt_estimate <- function(df, method, param = "gamma") {
  rownames(df) <- NULL
  structure(df, class = c("evt_estimate", "data.frame"), method = method,
            param = param)
}

#' Print an estimate path
#'
#' @param x An `evt_estimate` object.
#' @param digits Number of significant digits.
#' @param ... Passed to the data frame print method.
#' @return `x`, invisibly.
#' @export
print.evt_estimate <- function(x, digits = 4L, ...) {
  cat("<", attr(x, "method"), ">\n", sep = "")
  print(as.data.frame(x), digits = digits, row.names = FALSE, ...)
  invisible(x)
}

#' Plot an estimate path against k
#'
#' @param x An `evt_estimate` object with more than one value of k.
#' @param y Name of the column to plot; the main parameter by default.
#' @param ... Passed to [graphics::plot()].
#' @return `x`, invisibly.
#' @export
plot.evt_estimate <- function(x, y = attr(x, "param"), ...) {
  if (!y %in% names(x)) stop("no column `", y, "` to plot.", call. = FALSE)
  graphics::plot(x$k, x[[y]], type = "l", xlab = "k", ylab = y,
                 main = attr(x, "method"), ...)
  invisible(x)
}

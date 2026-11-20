test_that("the residual CV (equation 1) and T_m (equation 10)", {
  set.seed(1)
  x <- stats::rexp(64)
  z <- sort(x, decreasing = TRUE)
  cvk <- function(j) {
    top <- z[1:j]
    thr <- if (j == 64) 0 else z[j + 1]
    sqrt(mean(top^2) - mean(top)^2) / mean(top - thr)
  }
  tm <- sum(vapply(c(64, 32, 16, 8, 4), function(j) j * (cvk(j) - 1)^2, 0))
  res <- castillo2014(x, m = 4, nsim = 500)
  expect_equal(res$statistic, tm)
  # the whole-sample cv is the sample cv with divisor n
  expect_equal(cvk(64), sqrt(mean((x - mean(x))^2)) / mean(x))
  # Section 4.3: T4 from the printed CVs of the 256 largest values, 13.37
  cvs <- c(1.112, 1.100, 1.168, 1.184, 1.611)
  expect_equal(round(sum(c(256, 128, 64, 32, 16) * (cvs - 1)^2), 1), 13.4)
})

test_that("Table 1: the 'Approximate' row (a + b chi^2_nu)", {
  pub <- rbind(c(3.84, 6.63), c(6.93, 11.65), c(9.67, 15.94), c(12.18, 19.69),
               c(14.49, 23.01))
  for (m in 0:4) {
    ap <- castillo2014_approx_q <- castillo_approx(m)
    q <- ap[["a"]] + ap[["b"]] * stats::qchisq(c(.95, .99), ap[["nu"]])
    # the 99% point of chi-squared(1) is 6.635, printed 6.63 (truncated)
    expect_equal(trunc(q * 100) / 100, pub[m + 1, ], tolerance = 0.002)
  }
})

test_that("Table 1: the 'Asymptotic' row (Theorem 2) by simulation", {
  skip_on_cran()
  pub <- rbind(c(6.97, 11.62), c(9.60, 15.98), c(12.04, 19.69),
               c(14.39, 22.88))
  set.seed(2)
  for (m in 1:4) {
    idx <- 0:m
    lam <- eigen(2^(-abs(outer(idx, idx, "-")) / 2))$values
    s <- colSums(lam * matrix(stats::rnorm(2e5 * (m + 1))^2, m + 1))
    q <- stats::quantile(s, c(.95, .99), names = FALSE)
    expect_lt(max(abs(q / pub[m, ] - 1)), 0.025)
  }
})

test_that("Table 1: finite-sample points of T0; T1-T4 come out larger", {
  skip_on_cran()
  q <- function(n, m) {
    stats::quantile(evtests:::sim_null(paste("castillo", n, m), 2e4,
      function(ns) vapply(seq_len(ns), function(i) {
        castillo_stat(sort(stats::rexp(n), decreasing = TRUE), m)
      }, 0)), c(.95, .99), names = FALSE)
  }
  expect_lt(max(abs(q(100, 0) / c(3.36, 6.43) - 1)), 0.05)
  expect_lt(max(abs(q(500, 0) / c(3.71, 6.68) - 1)), 0.05)
  # T4 at n = 500: the 95% point (12.24) is reproduced; the 99% point
  # (22.86) is not (about 21.5 here), and the printed 99% points are not
  # monotone in n
  q4 <- q(500, 4)
  expect_lt(abs(q4[1] - 12.24), 0.25)
  expect_gt(abs(q4[2] - 22.86), 0.5)
})

test_that("Table 3: power of T0 against a Pareto tail, xi = 0.10", {
  skip_on_cran()
  set.seed(3)
  p <- replicate(1500, {
    x <- (stats::runif(500)^(-0.1) - 1) / 0.1
    castillo2014(x, m = 0, nsim = 1e4)$p_value
  })
  expect_lt(abs(mean(p < 0.05) - 0.610), 0.04)
})

test_that("the CV-plot", {
  set.seed(4)
  x <- stats::rexp(300)
  cv <- castillo2014_cv(x, k = c(10, 100, 299), level = 0.9)
  expect_equal(cv$upper, 1 + stats::qnorm(0.95) / sqrt(c(10, 100, 299)))
  expect_equal(cv$p_value, 1 - stats::pnorm(sqrt(cv$k) * (cv$cv - 1)))
  expect_equal(cv$cv[2], residual_cv(sort(x, decreasing = TRUE), 100))
})

test_that("validation", {
  expect_error(castillo2014(c(-1, 2, 3, 4)), "non-negative")
  expect_error(castillo2014(stats::rexp(20), m = 4), "n / 2\\^m")
  expect_error(castillo2014_cv(c(-1, 2, 3)), "non-negative")
})

test_that("the limit quantiles reproduce Table 1 of Drees, de Haan & Li", {
  tab <- evtests:::limit_tables$T
  ps <- c(0.995, 0.99, 0.975, 0.95, 0.9, 0.5, 0.1, 0.05, 0.025, 0.01, 0.005)
  gam <- c(2, 1.5, 1, 0.5, 0.25, 0, -0.25, -0.375, -0.45, -0.49)
  pub <- matrix(c(
    0.545, 0.513, 0.507, 0.525, 0.553, 0.621, 0.672, 0.739, 0.739, 0.889,
    0.477, 0.462, 0.459, 0.474, 0.494, 0.554, 0.604, 0.667, 0.726, 0.774,
    0.408, 0.389, 0.383, 0.390, 0.409, 0.459, 0.510, 0.558, 0.590, 0.641,
    0.349, 0.337, 0.330, 0.337, 0.355, 0.390, 0.431, 0.468, 0.500, 0.539,
    0.289, 0.281, 0.278, 0.285, 0.295, 0.318, 0.355, 0.381, 0.405, 0.435,
    0.151, 0.148, 0.147, 0.149, 0.154, 0.162, 0.178, 0.189, 0.199, 0.207,
    0.083, 0.082, 0.081, 0.082, 0.085, 0.089, 0.095, 0.099, 0.103, 0.105,
    0.071, 0.070, 0.070, 0.071, 0.073, 0.078, 0.080, 0.083, 0.087, 0.089,
    0.062, 0.062, 0.062, 0.063, 0.064, 0.068, 0.071, 0.073, 0.076, 0.077,
    0.053, 0.054, 0.054, 0.055, 0.056, 0.059, 0.060, 0.062, 0.066, 0.066,
    0.048, 0.049, 0.049, 0.050, 0.051, 0.052, 0.054, 0.055, 0.059, 0.058),
    nrow = length(ps), byrow = TRUE)
  sim <- vapply(gam, function(g) {
    stats::approx(evtests:::limit_tables$probs,
                  tab$q[, which(abs(tab$gamma - g) < 1e-9)], xout = ps)$y
  }, numeric(length(ps)))
  rel <- abs(sim / pub - 1)
  # the cell gamma = -0.45, p = .995 is printed as 0.739, a repetition of
  # the cell to its left; the simulated quantile is about 0.81
  expect_gt(sim[1, 9], 0.79)
  rel[1, 9] <- 0
  # the extreme quantiles (p = .995 and .005) carry more Monte Carlo error
  expect_lt(max(rel[-c(1, 11), ]), 0.035)
  expect_lt(max(rel[c(1, 11), ]), 0.06)
})

test_that("k T is the exact integral of (G - x)^2 / x", {
  set.seed(1)
  x <- stats::rnorm(400)
  k <- 40
  z <- sort(x, decreasing = TRUE)
  y <- z[1:k] - z[k + 1]
  fit <- grimshaw1993(y)
  v <- (1 + fit$gamma * y / fit$scale)^(-1 / fit$gamma)
  f <- function(u) vapply(u, function(s) (mean(v < s) - s)^2 / s, 0)
  brk <- sort(c(0, v, 1))
  num <- sum(vapply(seq_len(length(brk) - 1), function(j) {
    stats::integrate(f, brk[j], brk[j + 1], rel.tol = 1e-12)$value
  }, numeric(1)))
  res <- drees2006(x, k = k)
  expect_equal(res$statistic, k * num, tolerance = 1e-8)
  expect_equal(res$gamma, fit$gamma)
  # location and scale invariance
  expect_equal(drees2006(3 + 2 * x, k = k)$statistic, res$statistic,
               tolerance = 1e-8)
  # the exponential case gamma = 0 uses exp(-y / a)
  expect_equal(drees2006_stat(y, 0, 1),
               drees2006_stat(y, 1e-13, 1), tolerance = 1e-6)
})

test_that("Table 7 of Drees, de Haan & Li: empirical size, n = 1000", {
  skip_on_cran()
  set.seed(2)
  size <- function(qs, k, nrep = 1000) {
    p <- replicate(nrep, {
      e <- cumsum(stats::rexp(k + 1))
      s <- e / (e[k + 1] + stats::rgamma(1, 1000 - k))
      drees2006(qs(s), k = k)$p_value
    })
    mean(p[!is.na(p)] < 0.05)
  }
  cauchy <- function(s) 1 / tan(pi * s)
  ev0 <- function(s) -log(-log1p(-s))
  # printed (percent): Cauchy k = 50: 4.9, k = 300: 11.5; EV(0) k = 100: 4.3
  check <- function(sim, pub) {
    expect_lt(abs(100 * sim - pub), 3 * sqrt(pub * (100 - pub) / 1000) + 1)
  }
  check(size(cauchy, 50), 4.9)
  check(size(cauchy, 300), 11.5)
  check(size(ev0, 100), 4.3)
})

test_that("NA when there is no estimate or it is below -0.499", {
  expect_true(is.na(drees2006(c(4, 3, 2, 1), k = 3)$p_value))
  set.seed(3)
  r <- drees2006(stats::runif(2000), k = 200)    # gamma = -1
  expect_true(is.na(r$p_value) || r$gamma >= -0.499)
  tied <- drees2006(c(rep(5, 4), 1:3), k = 3)
  expect_true(is.na(tied$statistic))
  expect_s3_class(drees2006(stats::rexp(30)), "evt_test")
})

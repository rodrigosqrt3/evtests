test_that("generalized Pareto percentiles give back the shape and scale", {
  # excesses over Z[4M] are GPD(c, a); place Z[M], Z[2M] at the fitted 1/4 and
  # 1/2 upper percentiles of the excesses (equations 3.1-3.3)
  for (cc in c(-0.4, 0, 0.3, 1)) {
    a <- 2
    q <- function(p) if (cc == 0) -a * log(p) else a * (p^(-cc) - 1) / cc
    m <- 5
    z <- c(rep(10 + q(0.25), m), rep(10 + q(0.5), m), rep(10, 2 * m))
    z <- z + c(seq(0, 1e-9, length.out = m), rep(0, 3 * m))
    p <- pickands1975(z, k = m)
    expect_equal(p$gamma, cc, tolerance = 1e-6)
    expect_equal(p$scale, a, tolerance = 1e-6)
  }
})

test_that("equations 2.3 and 2.4 in closed form", {
  set.seed(2)
  x <- stats::rexp(100)
  z <- sort(x, decreasing = TRUE)
  m <- 12
  g <- log((z[m] - z[2 * m]) / (z[2 * m] - z[4 * m])) / log(2)
  a <- (z[2 * m] - z[4 * m]) /
    stats::integrate(function(s) exp(g * s), 0, log(2))$value
  p <- pickands1975(x, k = m)
  expect_equal(p$gamma, g)
  expect_equal(p$scale, a)
})

test_that("location and scale equivariance", {
  set.seed(3)
  x <- stats::rnorm(200)
  p <- pickands1975(x, k = c(5, 20, 50))
  q <- pickands1975(3 + 7 * x, k = c(5, 20, 50))
  expect_equal(q$gamma, p$gamma)
  expect_equal(q$scale, 7 * p$scale)
})

test_that("it is consistent for the three domains", {
  set.seed(4)
  n <- 40000
  m <- 1000
  # the asymptotic standard deviation is about 0.06 at M = 1000
  expect_lt(abs(pickands1975(stats::rexp(n), k = m)$gamma), 0.2)
  expect_lt(abs(pickands1975(stats::runif(n)^(-1 / 2), k = m)$gamma - 0.5),
            0.2)
  expect_lt(abs(pickands1975(1 - stats::runif(n)^(1 / 2), k = m)$gamma + 0.5),
            0.2)
})

test_that("structure, whole path and validation", {
  p <- pickands1975(1:40)
  expect_s3_class(p, "evt_estimate")
  expect_equal(p$k, 1:10)
  # equally spaced values: Z[M] - Z[2M] = M, Z[2M] - Z[4M] = 2M, so c = -1
  expect_equal(p$gamma, rep(-1, 10))
  expect_error(pickands1975(1:3), "at least 4")
  expect_error(pickands1975(1:40, k = 11), "between 1 and 10")
})

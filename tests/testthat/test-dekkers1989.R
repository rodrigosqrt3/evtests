test_that("equations 1.3, 1.5 and 1.7 in closed form", {
  set.seed(1)
  x <- stats::rlnorm(300)
  z <- sort(x, decreasing = TRUE)
  k <- 40
  d <- log(z[1:k]) - log(z[k + 1])
  m1 <- mean(d)
  m2 <- mean(d^2)
  res <- dekkers1989(x, k = k)
  expect_equal(res$M1, m1)
  expect_equal(res$M2, m2)
  expect_equal(res$gamma, m1 + 1 - 0.5 / (1 - m1^2 / m2))
  # M1 is the Hill estimator of gamma
  expect_equal(res$M1, hill1975(x, k = k)$gamma)
})

test_that("the variance of equation 3.10 and its special values", {
  expect_equal(moment_variance(0), 1)
  expect_equal(moment_variance(0.5), 1.25)
  # the gamma < 0 branch simplifies to
  # (1 - g)^2 (1 - 2g) (1 - g + 6 g^2) / ((1 - 3g)(1 - 4g))
  for (g in c(-0.1, -0.5, -1, -2)) {
    expect_equal(moment_variance(g),
                 (1 - g)^2 * (1 - 2 * g) * (1 - g + 6 * g^2) /
                   ((1 - 3 * g) * (1 - 4 * g)))
  }
  # continuous at 0
  expect_equal(moment_variance(-1e-8), 1, tolerance = 1e-6)
  # and it follows from the delta method on the matrix of equation 3.7 once
  # s22 has (1 - 3g)(1 - 4g) in the denominator
  g <- -0.3
  r1 <- 1 / (1 - g)
  r2 <- 2 / ((1 - g) * (1 - 2 * g))
  s <- matrix(c(1 / ((1 - g)^2 * (1 - 2 * g)),
                4 / ((1 - g)^2 * (1 - 2 * g) * (1 - 3 * g)),
                4 / ((1 - g)^2 * (1 - 2 * g) * (1 - 3 * g)),
                4 * (5 - 11 * g) / ((1 - g)^2 * (1 - 2 * g)^2 * (1 - 3 * g) *
                                      (1 - 4 * g))), 2)
  # gamma = 1 - 1 / (2 (1 - m1^2 / m2)) for gamma < 0 (M1 -> 0 at the scale
  # of f); gradient with respect to (m1, m2) at (r1, r2)
  q <- r1^2 / r2
  grad <- -0.5 / (1 - q)^2 * c(2 * r1 / r2, -r1^2 / r2^2)
  expect_equal(drop(t(grad) %*% s %*% grad), moment_variance(g))
})

test_that("the estimator has the variance of Corollary 3.2", {
  skip_on_cran()
  # generalized Pareto tails, top order statistics by the Renyi
  # representation, n large so that the bias of the logarithm is negligible
  set.seed(2)
  top <- function(n, k) {
    e <- cumsum(stats::rexp(k + 1))
    e / (e[k + 1] + stats::rgamma(1, n - k))
  }
  k <- 300
  sim <- function(q) replicate(2000, {
    dekkers1989(q(top(1e7, k)), k = k)$gamma
  })
  neg <- sim(function(p) 3 - p^0.5)          # gamma = -0.5
  pos <- sim(function(p) p^(-0.5))           # gamma = 0.5
  expect_equal(stats::var(sqrt(k) * neg), moment_variance(-0.5),
               tolerance = 0.1)
  expect_equal(stats::var(sqrt(k) * pos), moment_variance(0.5),
               tolerance = 0.1)
  expect_lt(abs(mean(neg) + 0.5), 0.02)
  expect_lt(abs(mean(pos) - 0.5), 0.02)
})

test_that("the end-point of equation 4.8", {
  set.seed(3)
  x <- 3 - stats::runif(5000)^0.5
  res <- dekkers1989(x, k = c(100, 300))
  z <- sort(x, decreasing = TRUE)
  for (i in 1:2) {
    kk <- res$k[i]
    expect_equal(res$endpoint[i],
                 z[kk + 1] * res$M1[i] * (1 - 1 / res$gamma[i]) + z[kk + 1])
  }
  expect_equal(res$endpoint, c(3, 3), tolerance = 0.02)
  pareto <- dekkers1989(stats::runif(500)^(-1), k = 50)
  expect_gt(pareto$gamma, 0)
  expect_true(is.na(pareto$endpoint))
})

test_that("scale invariance, but not location invariance", {
  set.seed(4)
  x <- stats::rexp(300) + 5
  expect_equal(dekkers1989(10 * x, k = 30)$gamma, dekkers1989(x, k = 30)$gamma)
  expect_false(isTRUE(all.equal(dekkers1989(x + 100, k = 30)$gamma,
                                dekkers1989(x, k = 30)$gamma)))
})

test_that("structure and validation", {
  res <- dekkers1989(us_cities_1940$population)
  expect_s3_class(res, "evt_estimate")
  expect_equal(res$k, 2:29)
  expect_error(dekkers1989(1:2), "at least 3")
  expect_error(dekkers1989(c(5, 4, 3, -1), k = 3), "positive")
  expect_error(dekkers1989(1:10, k = 1), "between 2 and 9")
  expect_true(is.na(moment_variance(NaN)))
})

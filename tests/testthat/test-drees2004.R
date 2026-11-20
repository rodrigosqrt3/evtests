test_that("it is the Grimshaw fit to the excesses over X[n-k,n]", {
  set.seed(1)
  x <- stats::rnorm(300)
  d <- drees2004(x, k = c(20, 60))
  z <- sort(x, decreasing = TRUE)
  for (i in 1:2) {
    kk <- d$k[i]
    fit <- grimshaw1993(z[1:kk] - z[kk + 1])
    expect_equal(d$gamma[i], fit$gamma)
    expect_equal(d$scale[i], fit$scale)
  }
})

test_that("standard errors follow Corollary 2.1", {
  set.seed(2)
  d <- drees2004(stats::runif(2000)^(-0.25), k = 200)
  g <- d$gamma
  expect_equal(d$se_gamma, (1 + g) / sqrt(200))
  expect_equal(d$se_scale_rel, sqrt((2 + 2 * g + g^2) / 200))
})

test_that("the estimator is close to normal with the variance of Sigma", {
  skip_on_cran()
  # GPD sample with gamma = 0.25: no bias term, variance (1 + gamma)^2 / k
  set.seed(3)
  k <- 200
  g <- replicate(500, drees2004((stats::runif(1000)^(-0.25) - 1) / 0.25,
                                k = k)$gamma)
  expect_lt(abs(mean(g) - 0.25), 3 * 1.25 / sqrt(k) / sqrt(500) + 0.01)
  expect_equal(stats::sd(g), 1.25 / sqrt(k), tolerance = 0.12)
})

test_that("gamma <= -1/2 and samples without an estimate", {
  # excesses whose estimate is the boundary point gamma = -1 (see the tests
  # of grimshaw1993): the asymptotic result does not apply
  u <- c(0.659, 0.199, 0.05, 0.347, 0.249, 0.655, 0.273, 0.034, 0.343, 0.155,
         0.211, 0.163, 0.002)
  d <- drees2004(c(1 + u, 1), k = 13)
  expect_equal(d$gamma, -1)
  expect_true(is.na(d$se_gamma) && is.na(d$se_scale_rel))
  # tied excesses: all equal to zero
  tied <- drees2004(c(rep(5, 4), 1:3), k = 3)
  expect_false(tied$exists)
  expect_true(is.na(tied$gamma))
  # no stationary maximum
  none <- drees2004(c(4, 3, 2, 1), k = 3)
  expect_false(none$exists)
})

test_that("structure and validation", {
  d <- drees2004(us_cities_1940$population)
  expect_s3_class(d, "evt_estimate")
  expect_equal(d$k, 2:29)
  expect_error(drees2004(1:2), "at least 3")
  expect_error(drees2004(1:10, k = 1), "between 2 and 9")
})

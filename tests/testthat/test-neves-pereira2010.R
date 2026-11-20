test_that("equations 5, 6, 17 and 18 in closed form", {
  set.seed(1)
  x <- stats::rexp(500) + 2
  z <- sort(x, decreasing = TRUE)
  n <- 500
  k <- 40
  d <- log(z[1:k]) - log(z[k + 1])
  a <- z[k + 1] * (mean(d) / 2) / (1 - mean(d)^2 / mean(d^2))
  t1 <- mean(z[2:(k + 1)] - z[k + 1] - a) / (z[1] - z[k + 1])
  r1 <- neves_pereira2010(x, k = k)
  expect_equal(r1$statistic, sqrt(k) * log(k) * t1)
  expect_equal(r1$p_value, 2 * stats::pnorm(-abs(r1$statistic)))
  t2 <- sum((1:(k - 1)) * (z[1:(k - 1)] - z[2:k])) / (k * z[k + 1])
  r2 <- neves_pereira2010(x, k = k, statistic = "T2")
  expect_equal(r2$statistic, sqrt(k) * (log(n / k) * t2 - 1))
})

test_that("Figure 1: rejection rates at k = 100, n = 1000", {
  skip_on_cran()
  # read from the figure: exponential about .04 (T1*) and .08 (T2*);
  # GPD(-0.1) about .25 and .55; GPD(-0.3) about .97 and 1
  rate <- function(gen, nrep = 1500) {
    rowMeans(replicate(nrep, {
      x <- gen()
      c(neves_pereira2010(x, k = 100)$p_value < 0.05,
        neves_pereira2010(x, k = 100, statistic = "T2")$p_value < 0.05)
    }))
  }
  # standard models, without shifts: the statistics are not location
  # invariant
  gpd <- function(g) function() ((stats::runif(1000))^(-g) - 1) / g
  set.seed(2)
  ex <- rate(function() stats::rexp(1000))
  expect_lt(max(abs(ex - c(.04, .08))), 0.04)
  m1 <- rate(gpd(-0.1))
  expect_lt(max(abs(m1 - c(.25, .55))), 0.06)
  m3 <- rate(gpd(-0.3))
  expect_gt(m3[1], 0.88)
  expect_gt(m3[2], 0.98)
})

test_that("T2* shifts to the left for the negative Frechet example", {
  skip_on_cran()
  # negative Frechet with x0 = beta = 1: U(t) = 1 - 1 / log(t), Gumbel
  # domain with end-point 1; the paper says T2* shifts to the right
  set.seed(3)
  t2 <- replicate(300, {
    x <- 1 + 1 / log1p(-stats::runif(1000))
    neves_pereira2010(x, k = 100, statistic = "T2")$statistic
  })
  expect_lt(mean(t2), -3)
})

test_that("validation", {
  expect_error(neves_pereira2010(c(5, 4, 3, -1), k = 3), "positive")
  expect_error(neves_pereira2010(c(1, 1, 1, 0.5), k = 2), "not all be equal")
  expect_s3_class(neves_pereira2010(stats::rexp(30) + 1), "evt_test")
})

test_that("equations 6 and 7 of Fraga Alves et al. (2016) agree", {
  set.seed(1)
  x <- stats::rnorm(200)
  z <- sort(x, decreasing = TRUE)
  k <- 30
  i <- 0:(k - 1)
  eq6 <- z[1] + z[k + 1] - sum(log(1 + 1 / (k + i)) * z[k + i + 1]) / log(2)
  a <- log((k + i + 1) / (k + i)) / log(2)
  eq7 <- z[1] + sum(a * (z[k + 1] - z[k + i + 1]))
  expect_equal(sum(a), 1)
  expect_equal(eq6, eq7)
  res <- fraga_alves_neves2014(x, k = k)
  expect_equal(res$endpoint, eq7)
  expect_gte(res$endpoint, max(x))
  # location and scale equivariance
  expect_equal(fraga_alves_neves2014(2 + 3 * x, k = k)$endpoint,
               2 + 3 * res$endpoint)
})

test_that("Table 1 of Fraga Alves et al. (2016): MAX and FAN, model 4", {
  skip_on_cran()
  # Beta(1, -1/gamma) with end-point 1, n = 1000, 300 replicates; average
  # absolute error at the best k (k* = 2k top order statistics for FAN)
  l1 <- function(g) {
    ks <- c(seq(2, 50, by = 2), seq(60, 250, by = 10))
    err <- replicate(300, {
      x <- 1 - stats::runif(1000)^(-g)
      c(max = 1 - max(x),
        fan = abs(fraga_alves_neves2014(x, k = ks)$endpoint - 1))
    })
    c(max = mean(err[1, ]), fan = min(rowMeans(err[-1, ])))
  }
  set.seed(2)
  half <- l1(-1 / 2)
  fifth <- l1(-1 / 5)
  # printed: MAX .029 and .234, FAN .014 and .171
  expect_equal(unname(half), c(.029, .014), tolerance = 0.15)
  expect_equal(unname(fifth), c(.234, .171), tolerance = 0.1)
})

test_that("the upper bound is equation 14 with the estimators 19-21", {
  set.seed(4)  # a sample whose estimated index (-0.41) is in (-1/2, 0)
  x <- 1 - stats::runif(2000)^0.4
  k <- 100
  res <- fraga_alves_neves2014(x, k = k, bound = TRUE, alpha = 0.1)
  z <- sort(x, decreasing = TRUE)
  e <- z[1:k] - z[k + 1]
  n1 <- mean(e)
  n2 <- mean(e^2)
  g <- 1 - 1 / (2 * (1 - n1^2 / n2))
  a <- n1 * (1 - g)
  h <- (1 / g) * ((2^(-g) - 1) / (g * log(2)) + 1)
  q <- (-log(0.1))^(-g) / g
  expect_equal(res$gamma, g)
  expect_equal(res$scale, a)
  expect_equal(res$upper, res$endpoint - a * (h + k^g * q))
  expect_gt(res$upper, res$endpoint)
})

test_that("coverage of the bound for gamma = -0.4", {
  skip_on_cran()
  set.seed(4)
  hit <- replicate(300, {
    x <- 1 - stats::runif(1000)^0.4
    r <- suppressWarnings(fraga_alves_neves2014(x, k = 100, bound = TRUE))
    r$upper >= 1
  })
  expect_gt(mean(hit, na.rm = TRUE), 0.8)
})

test_that("the bound is NA outside (-1/2, 0), with a warning", {
  set.seed(5)
  expect_warning(r <- fraga_alves_neves2014(stats::runif(200), k = c(10, 20),
                                            bound = TRUE),
                 "not in \\(-1/2, 0\\)")
  expect_true(all(is.na(r$upper)))
  expect_error(fraga_alves_neves2014(1:10, k = 1, bound = TRUE), "between 2")
  expect_error(fraga_alves_neves2014(1:10, bound = TRUE, alpha = 2), "alpha")
  expect_equal(fraga_alves_neves2014(1:10)$k, 1:5)
})

test_that("the interval of Theorem 3 of Fraga Alves & Neves (2014)", {
  set.seed(8)
  x <- 1 - (-log(stats::runif(1000)))^(-1 / 0.5) # negative Frechet, beta 0.5
  r <- fraga_alves_neves2014(x, k = 48, interval = TRUE, alpha = 0.1)
  z <- sort(x, decreasing = TRUE)
  ah <- mean(z[1:48] - z[49])
  expect_equal(r$mean_excess, ah)
  g <- -log(-log(c(0.95, 0.05))) - log(2) / 2
  expect_equal(c(r$lower_gumbel, r$upper_gumbel), r$endpoint - ah * g)
  expect_lt(r$lower_gumbel, r$endpoint)
})

test_that("coverage of the Theorem 3 interval, negative Frechet", {
  skip_on_cran()
  # Section 5: negative Frechet with end-point 1, n = 1000; beta < 1
  # satisfies the second-order condition; k = (log n)^2 = 48 (p. 1821)
  set.seed(9)
  hit <- replicate(500, {
    x <- 1 - (-log(stats::runif(1000)))^(-1 / 0.5)
    r <- fraga_alves_neves2014(x, k = 48, interval = TRUE)
    r$lower_gumbel <= 1 && 1 <= r$upper_gumbel
  })
  # conservative: all 500 intervals cover the end-point
  expect_gt(mean(hit), 0.95)
})

test_that("G* is equation 17 and has a Gumbel limit for exponential tails", {
  set.seed(6)
  x <- stats::rexp(500)
  z <- sort(x, decreasing = TRUE)
  k <- 40
  g <- (fraga_alves_neves2014(x, k = k)$endpoint - z[k + 1]) /
    (z[k + 1] - z[2 * k + 1])
  res <- fraga_alves2016(x, k = k, alternative = "positive", nsim = 2000)
  expect_equal(res$statistic, log(2) * g - (log(k) + log(2) / 2))
  expect_equal(res$p_value, 1 - exp(-exp(-res$statistic)))
  # the simulated exponential null is close to Gumbel for large k
  nul <- sim_null(paste("fanr", 200), 4000, function(ns) {
    vapply(seq_len(ns), function(i) {
      fanr_stat(c(sort(stats::rexp(400), decreasing = TRUE), 0), 200)
    }, numeric(1))
  })
  q <- stats::quantile(nul, c(.1, .5, .9), names = FALSE)
  expect_lt(max(abs(q - (-log(-log(c(.1, .5, .9)))))), 0.25)
})

test_that("G* rejects in the direction of the sign of gamma", {
  set.seed(7)
  expect_lt(fraga_alves2016(stats::runif(2000), k = 200,
                            alternative = "negative", nsim = 500)$p_value,
            0.01)
  expect_lt(fraga_alves2016(stats::runif(2000)^(-1), k = 200,
                            alternative = "positive", nsim = 500)$p_value,
            0.01)
  expect_error(fraga_alves2016(c(3, 2, 1, 1, 1), k = 2), "must differ")
  expect_equal(fraga_alves2016(stats::rexp(11), nsim = 200)$k, 1:5)
})

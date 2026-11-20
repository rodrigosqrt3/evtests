test_that("equations 1.3, 3.1 and 3.2 with Hasofer and Wang's W", {
  set.seed(1)
  x <- stats::rnorm(400)
  z <- sort(x, decreasing = TRUE)
  k <- 60
  e <- z[1:k] - z[k + 1]
  r <- mean(e^2) / mean(e)^2
  rs <- neves_fraga_alves2007(x, k = k, alternative = "positive")
  expect_equal(rs$statistic, sqrt(k / 4) * (r - 2))
  expect_equal(rs$p_value, 1 - stats::pnorm(rs$statistic))
  ws <- neves_fraga_alves2007(x, k = k, statistic = "W",
                              alternative = "negative")
  w <- hasofer_wang1992(x, k = k, null = "normal")$statistic
  expect_equal(ws$statistic, sqrt(k / 4) * (k * w - 1))
  expect_equal(ws$p_value, 1 - stats::pnorm(ws$statistic))
  two <- neves_fraga_alves2007(x, k = k)
  expect_equal(two$p_value, 2 * min(stats::pnorm(rs$statistic),
                                    1 - stats::pnorm(rs$statistic)))
})

test_that("the limits of Theorem 3.1 and Remark 3.1", {
  # c2 / c1^2 = 2 + 2g / (1 - 2g) and c1^2 / (c2 - c1^2) = 1 - 2g
  for (g in c(-0.5, -0.1, 0.2)) {
    c1 <- 1 / (1 - g)
    c2 <- 2 / ((1 - g) * (1 - 2 * g))
    expect_equal(c2 / c1^2, 2 + 2 * g / (1 - 2 * g))
    expect_equal(c1^2 / (c2 - c1^2), 1 - 2 * g)
  }
})

test_that("Table 1 of Neves & Fraga Alves (2007) is reproduced", {
  skip_on_cran()
  # generalized Pareto, n = 1000, two-sided tests at 5%
  rate <- function(g, k, nrep = 5000) {
    rej <- replicate(nrep, {
      top <- top_gpd(1000, k, g)
      c(neves_fraga_alves2007(top, k = k)$p_value < 0.05,
        neves_fraga_alves2007(top, k = k, statistic = "W")$p_value < 0.05)
    })
    rowMeans(rej)
  }
  set.seed(2)
  # columns R* and W*; both the published and the simulated rates have a
  # Monte Carlo error (5000 replicates each), so their difference has
  # standard error sqrt(2) se
  check <- function(sim, pub) {
    se <- sqrt(2 * pub * (1 - pub) / 5000)
    expect_true(all(abs(sim - pub) < 3 * se + 0.003))
  }
  check(rate(0, 50), c(.026, .051))
  check(rate(0, 100), c(.034, .051))
  check(rate(-0.2, 50), c(.020, .283))
  check(rate(-0.2, 100), c(.165, .508))
  check(rate(0.2, 100), c(.524, .323))
  check(rate(0.4, 50), c(.655, .384))
})

test_that("the printed W (equation 1.2) does not reproduce the table", {
  skip_on_cran()
  set.seed(3)
  k <- 50
  rej <- replicate(10000, {
    e <- top_gpd(1000, k, 0)
    e <- e[1:k] - e[k + 1]
    w <- (1 / k) * mean(e)^2 / (mean(e^2) - mean(e)^2)
    abs(sqrt(k / 4) * (k * w - 1)) > stats::qnorm(0.975)
  })
  # published size .051; the printed formula gives about .072
  expect_gt(mean(rej), 0.062)
})

test_that("validation", {
  expect_error(neves_fraga_alves2007(c(1, 1, 1, 0), k = 2), "not all be equal")
  expect_error(neves_fraga_alves2007(1:10, statistic = "Z"), "should be one")
  expect_equal(neves_fraga_alves2007(1:10)$k, 2:9)
})

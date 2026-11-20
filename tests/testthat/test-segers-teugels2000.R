test_that("the score has mean 0 and variance 1/5 under the null law", {
  # Corollary 2.2: P(xi <= x) = 1 - 2 / (1 + x), x > 1
  set.seed(1)
  xi <- 2 / (1 - stats::runif(2e5)) - 1
  expect_lt(abs(mean(st_score(xi))), 0.003)
  expect_equal(stats::var(st_score(xi)), 1 / 5, tolerance = 0.01)
  # Galton's ratio of an exponential sample has this law
  g <- replicate(2e4, {
    b <- sort(stats::rexp(20), decreasing = TRUE)[1:3]
    (b[1] - b[3]) / (b[2] - b[3])
  })
  expect_lt(abs(mean(g <= 3) - (1 - 2 / 4)), 0.01)
})

test_that("S and Z in closed form", {
  set.seed(2)
  x <- stats::rnorm(200)
  xi <- vapply(1:10, function(i) {
    b <- sort(x[(i - 1) * 20 + 1:20], decreasing = TRUE)[1:3]
    (b[1] - b[3]) / (b[2] - b[3])
  }, 0)
  tt <- 1 - 6 * xi / (1 + xi)^2
  s <- segers_teugels2000(x)
  expect_equal(s$statistic, (5 / 10) * sum(tt)^2)
  expect_equal(s$p_value, 1 - stats::pchisq(s$statistic, 1))
  z <- segers_teugels2000(x, alternative = "negative")
  expect_equal(z$statistic, sqrt(5 / 10) * sum(tt))
  expect_equal(z$p_value, stats::pnorm(z$statistic))
  expect_equal(segers_teugels2000(3 + 2 * x)$statistic, s$statistic)
  zp <- segers_teugels2000(x, alternative = "positive")
  expect_equal(zp$p_value, 1 - stats::pnorm(zp$statistic))
})

test_that("Section 4: size for GPD(0) and the direction of the alternatives", {
  skip_on_cran()
  set.seed(3)
  p <- replicate(1000, segers_teugels2000(stats::rexp(1000))$p_value)
  # the paper's figures show sizes close to nominal (our runs: .111, .052,
  # .006 at .1, .05, .01)
  expect_lt(abs(mean(p < 0.05) - 0.05), 0.02)
  zp <- replicate(200, segers_teugels2000(stats::runif(1000)^(-0.5),
                                          alternative = "positive")$p_value)
  expect_gt(mean(zp < 0.05), 0.3)
  zn <- replicate(200, segers_teugels2000(1 - stats::runif(1000)^0.5,
                                          alternative = "negative")$p_value)
  expect_gt(mean(zn < 0.05), 0.3)
})

test_that("validation", {
  expect_error(segers_teugels2000(1:30), "two subsamples")
  expect_error(segers_teugels2000(1:100, block_size = 2), "at least 3")
  expect_error(segers_teugels2000(c(rep(c(3, 1, 1), 20))), "tied")
  expect_error(segers_teugels2000(c(1, NA)), "numeric")
})

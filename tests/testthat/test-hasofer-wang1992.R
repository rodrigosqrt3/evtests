test_that("W, its excess form and the Greenwood relation of Section 7", {
  set.seed(1)
  x <- stats::rnorm(200)
  z <- sort(x, decreasing = TRUE)
  k <- 30
  top <- z[1:(k + 1)]
  w <- (k + 1) * (mean(top) - top[k + 1])^2 / (k * sum((top - mean(top))^2))
  res <- hasofer_wang1992(x, k = k, null = "normal")
  expect_equal(res$statistic, w)
  # in terms of the k excesses: W = 1 / ((k + 1) R - k), R = M2 / M1^2
  e <- top[1:k] - top[k + 1]
  r <- mean(e^2) / mean(e)^2
  expect_equal(w, 1 / ((k + 1) * r - k))
  # Greenwood: script-W = 1/K + 1/(K (K - 1) W) = sum(Q^2), K = k + 1
  expect_equal(1 / (k + 1) + 1 / ((k + 1) * k * w), sum((e / sum(e))^2))
  # location and scale invariance
  expect_equal(hasofer_wang1992(5 + 3 * x, k = k, null = "normal")$statistic,
               w)
})

test_that("Table 7: the theoretical mean and standard deviation", {
  # the paper's k counts the k + 1 top order statistics: k = 20, ..., 60
  kk <- c(20, 30, 40, 50, 60) - 1
  # 1 / 59 = .016949 is printed as .0170
  expect_equal(round(1 / kk, 4), c(.0526, .0345, .0256, .0204, .0169))
  expect_equal(round(sqrt(4 * (kk - 1) / (kk^2 * (kk + 2) * (kk + 3))), 4),
               c(.0208, .0116, .0076, .0055, .0042))
})

test_that("Table 8: percentage points of the limiting null", {
  skip_on_cran()
  tab <- rbind(
    `13` = c(333.8, 395.6, 457.7, 538.4, 1599.5, 1827.5, 1998.5, 2204.6),
    `30` =c(165.0, 189.2, 212.6, 241.7, 535.5, 593.6, 650.7, 719.9),
    `40` = c(131.7, 149.5, 166.1, 186.4, 377.1, 413.6, 451.5, 497.0),
    `50` = c(110.4, 124.1, 136.6, 151.8, 289.5, 316.5, 340.2, 371.1),
    `60` = c(96.4, 106.9, 116.6, 128.8, 232.9, 253.8, 270.8, 293.3),
    `100` = c(63.6, 69.2, 74.3, 80.2, 129.6, 138.4, 146.4, 155.5),
    `200` = c(35.4, 37.6, 39.9, 41.9, 59.4, 62.2, 64.8, 68.1)
  )
  probs <- c(.01, .025, .05, .10, .90, .95, .975, .99)
  dev <- function(K) {
    q <- stats::quantile(hw_null(as.integer(K) - 1L, 4e5), probs,
                         names = FALSE)
    1e4 * q / tab[K, ] - 1
  }
  # 30 to 100 order statistics: within 1.5% in every column
  for (K in c("30", "40", "50", "60", "100")) {
    expect_lt(max(abs(dev(K))), 0.015)
  }
  # 200: the table is 1-2% below the simulated points in every column, and
  # for 13 order statistics 2-3% off and 7% too small at the upper 1% point
  # (2204.6 against about 2370); the table was interpolated from Stephens
  # (1981)
  expect_lt(max(abs(dev("200"))), 0.025)
  d13 <- dev("13")
  expect_lt(max(abs(d13[1:7])), 0.035)
  expect_gt(d13[8], 0.05)
})

test_that("Tables 1 and 2: n = 300, 50 top order statistics, alpha = 0", {
  skip_on_cran()
  set.seed(2)
  w <- replicate(5000, {
    x <- -log(-log(stats::runif(300)))                   # Gumbel
    hw_stat(sort(x, decreasing = TRUE)[1:50])
  })
  expect_equal(round(mean(w), 4), .0228, tolerance = 0.0004 / .0228)
  expect_equal(round(stats::sd(w), 4), .0057, tolerance = 0.0004 / .0057)
  q <- 1e4 * stats::quantile(w, c(.05, .10, .90, .95), names = FALSE)
  expect_lt(max(abs(q / c(145, 159, 301, 328) - 1)), 0.03)
})

test_that("Tables 1 and 3: GEV alternatives", {
  skip_on_cran()
  set.seed(3)
  gev <- function(a) {
    x <- -log(-log(stats::runif(300)))
    (exp(a * x) - 1) / a
  }
  sim <- function(a) replicate(3000, {
    hw_stat(sort(gev(a), decreasing = TRUE)[1:50])
  })
  # Table 1: mean and standard deviation (x 1e4)
  tab1 <- rbind(c(1, 38, 30), c(1 / 2, 92, 48), c(1 / 4, 147, 55),
                c(-1 / 4, 320, 63), c(-1 / 2, 420, 80), c(-1, 620, 128))
  # Table 3: power of the one-sided test with the critical values of
  # Table 2 (lower 116, 145, 159, 180; x 1e-4), at levels .01 .05 .1 .2
  tab3 <- list(`1` = c(.975, .990, .995, .999), `0.5` = c(.720, .844, .895,
                                                          .947),
               `0.25` = c(.300, .500, .590, .720))
  for (i in seq_len(nrow(tab1))) {
    w <- sim(tab1[i, 1])
    se <- tab1[i, 3] / sqrt(3000)
    expect_lt(abs(1e4 * mean(w) - tab1[i, 2]), 4 * se + 1)
    expect_equal(1e4 * stats::sd(w), tab1[i, 3], tolerance = 0.06)
    key <- as.character(tab1[i, 1])
    if (key %in% names(tab3)) {
      pw <- vapply(c(116, 145, 159, 180) / 1e4, function(cv) mean(w < cv), 0)
      expect_lt(max(abs(pw - tab3[[key]])), 0.035)
    }
  }
  # for negative alpha the published powers are lower than those obtained
  # with the critical values of Table 2 (alpha = -1/2, level .01: .600
  # printed, about .69 here), although the distribution of W agrees with
  # Table 1; the critical values used for those rows are not stated
  w <- sim(-1 / 2)
  expect_gt(mean(w > .0378) - .600, 0.05)
})

test_that("p-values: simulated null, normal approximation, alternatives", {
  set.seed(4)
  x <- stats::rexp(500)
  k <- 40
  sim <- hasofer_wang1992(x, k = k, alternative = "positive")
  nul <- hw_null(k, 1e5)
  expect_equal(sim$p_value, mean(nul <= sim$statistic))
  neg <- hasofer_wang1992(x, k = k, alternative = "negative")
  expect_equal(neg$p_value, mean(nul >= sim$statistic))
  two <- hasofer_wang1992(x, k = k)
  expect_equal(two$p_value, min(1, 2 * min(sim$p_value, neg$p_value)))
  nor <- hasofer_wang1992(x, k = k, null = "normal", alternative = "negative")
  s <- sqrt(4 * (k - 1) / (k^2 * (k + 2) * (k + 3)))
  expect_equal(nor$p_value, 1 - stats::pnorm((sim$statistic - 1 / k) / s))
  # the simulated null does not change the user's random stream
  set.seed(5)
  a <- stats::runif(1)
  set.seed(5)
  null_cache$tmp <- NULL
  hasofer_wang1992(x, k = 23, nsim = 1000)
  expect_equal(stats::runif(1), a)
})

test_that("it rejects for uniform and Pareto tails, in the right direction", {
  set.seed(6)
  expect_lt(hasofer_wang1992(stats::runif(1000), k = 100,
                             alternative = "negative")$p_value, 0.01)
  expect_lt(hasofer_wang1992(stats::runif(1000)^(-1), k = 100,
                             alternative = "positive")$p_value, 0.01)
})

test_that("structure, printing, plotting and validation", {
  set.seed(7)
  res <- hasofer_wang1992(stats::rexp(60), nsim = 1000)
  expect_s3_class(res, "evt_test")
  expect_equal(res$k, 2:59)
  expect_output(print(res), "gamma != 0")
  expect_output(print(hasofer_wang1992(stats::rexp(60), k = 9,
                                       alternative = "positive")),
                "gamma > 0")
  expect_output(print(hasofer_wang1992(stats::rexp(60), k = 9,
                                       alternative = "negative")),
                "gamma < 0")
  pdf(NULL)
  on.exit(dev.off())
  expect_invisible(plot(res))
  expect_invisible(plot(res, "statistic"))
  expect_error(plot(res, "nothing"), "no column")
  expect_error(hasofer_wang1992(c(1, 1, 1, 0), k = 2), "not all be equal")
  expect_error(hasofer_wang1992(1:10, nsim = 10), "nsim")
  expect_error(hasofer_wang1992(1:10, alternative = "up"), "should be one")
})

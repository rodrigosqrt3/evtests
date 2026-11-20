test_that("Example 1: bushfires, Frechet domain (log scale), k = 27 top", {
  w <- wang1995(bushfires$area, alpha = 0.1, transform = "log", nsim = 4e4)
  expect_equal(attr(w, "selected"), 26L)      # 27 top observations
})

test_that("Example 2: mice, Gumbel domain at alpha = 20%, k = 8 top", {
  w <- wang1995(-mice_survival$weeks, alpha = 0.2, nsim = 4e4)
  expect_equal(attr(w, "selected"), 7L)       # 8 top observations
})

test_that("Example 1 on the original scale gives 13 top (12 printed)", {
  # W for 13 top observations (k = 12) is just inside the acceptance region
  # (lower 5% point 0.5816 by simulation, 0.5950 from Table 8 of Hasofer &
  # Wang, both times 13; W x 13 = 0.5956); the paper's 12 implies a
  # rejection there
  w <- wang1995(bushfires$area, alpha = 0.1, nsim = 4e4)
  expect_equal(attr(w, "selected"), 12L)      # 13 top observations
  # Wang's W is 13 times Hasofer and Wang's
  expect_equal(round(13 * w$statistic[w$k == 12], 4), 0.5956)
})

test_that("the rule: last k before the first rejection", {
  set.seed(1)
  x <- stats::rexp(200)
  w <- wang1995(x, alpha = 0.2, kmax = 60, nsim = 2000)
  rej <- which(w$p_value < 0.2)
  expect_equal(attr(w, "selected"),
               if (length(rej)) w$k[rej[1]] - 1L else 60L)
  expect_equal(w$k, 2:60)
  # the end-point transform
  u <- stats::runif(100)
  we <- wang1995(u, transform = "endpoint", q = 1, kmax = 20, nsim = 1000)
  expect_equal(we$statistic,
               hasofer_wang1992(-log(1 - u), k = 2:20, nsim = 1000)$statistic)
})

test_that("validation", {
  expect_error(wang1995(c(-1, 2, 3), transform = "log"), "positive")
  expect_error(wang1995(1:10, transform = "endpoint", q = 5), "larger")
  expect_error(wang1995(1:10, alpha = 0), "alpha")
})

tm <- olympic_100m$time

test_that("Section 3 of Hall & Wang (1999): the 100 m example", {
  res <- hall_wang1999(tm, k = 11, alpha = 0.1, side = "lower")
  # printed: theta = 9.91 (9.9158; the paper truncates), alpha = 0.82
  # (0.826), 90% interval (9.80, 9.92)
  expect_equal(trunc(res$endpoint * 100) / 100, 9.91)
  expect_equal(round(res$endpoint, 4), 9.9158)
  expect_equal(trunc(res$shape * 100) / 100, 0.82)
  expect_equal(round(res$lower, 2), 9.80)
  expect_equal(round(res$upper, 2), 9.92)
  expect_lte(res$upper, 9.92)
  # G at the estimate and at -Inf (notes: .16493 and .27195)
  expect_equal(round(res$statistic, 5), 0.16493)
})

test_that("G(theta) at the values of Figure 2", {
  g <- sort(tm) - min(tm)
  d <- function(theta) min(tm) - theta
  vals <- vapply(c(9.70, 9.75, 9.80, 9.85, 9.90), function(th) {
    hw_g(d(th), g)
  }, numeric(1))
  expect_equal(round(vals, 4), c(0.2312, 0.2247, 0.2156, 0.2017, 0.1762))
  expect_equal(round(hw_g_inf(g), 5), 0.27195)
})

test_that("the interval is where G is below the null quantile", {
  res <- hall_wang1999(tm, k = 11, alpha = 0.1, side = "lower")
  g <- sort(tm) - min(tm)
  expect_equal(hw_g(min(tm) - res$lower, g), res$critical, tolerance = 1e-6)
  nul <- sim_null(paste("greenwood", 11L), 1e5,
                  function(ns) greenwood_sim(11L, ns))
  expect_equal(res$critical, stats::quantile(nul, 0.9, names = FALSE))
})

test_that("upper and lower end-points are mirror images", {
  set.seed(1)
  x <- 3 - stats::rexp(30)^(1 / 3)        # Weibull tail, end-point 3
  up <- hall_wang1999(x, k = 10)
  lo <- hall_wang1999(-x, k = 10, side = "lower")
  expect_equal(up$endpoint, -lo$endpoint)
  expect_equal(c(up$lower, up$upper), -c(lo$upper, lo$lower))
  expect_gte(up$endpoint, max(x))
})

test_that("coverage under the exact limiting model", {
  skip_on_cran()
  set.seed(2)
  cover <- replicate(400, {
    x <- 5 - cumsum(stats::rexp(12))^(1 / 3)   # theta = 5, alpha = 3
    r <- suppressWarnings(hall_wang1999(x, k = 11, nsim = 2e4))
    isTRUE(r$lower <= 5 && 5 <= r$upper)    # an empty region does not cover
  })
  expect_lt(abs(mean(cover) - 0.95), 0.03)
})

test_that("the shape from the GLS slope", {
  set.seed(3)
  d <- cumsum(stats::rexp(200))^(1 / 2.5)    # alpha = 2.5
  expect_equal(hw_shape(d), 2.5, tolerance = 0.1)
})

test_that("validation and a record without an estimate", {
  expect_error(hall_wang1999(1:3), "at least 4")
  expect_error(hall_wang1999(c(1, 1, 1, 1, 0), k = 3), "not all be equal")
  r <- hall_wang1999(c(10, 9, 7, 4, 0), k = 4, nsim = 2e3)
  expect_s3_class(r, "evt_estimate")
})

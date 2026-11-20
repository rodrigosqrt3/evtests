y <- wheaton_river$excess

test_that("Table 3: estimates and statistics as the smallest are deleted", {
  stopifnot(length(y) == 72, abs(sum(y) - 878.7) < 1e-9)
  pub <- rbind(c(72, -0.006, 12.14, 0.2309, 1.452),
               c(71, 0.019, 12.63, 0.208, 1.317),
               c(70, 0.046, 13.15, 0.187, 1.197),
               c(69, 0.047, 13.34, 0.162, 1.039),
               c(68, 0.077, 13.93, 0.143, 0.967),
               c(67, 0.078, 14.13, 0.125, 0.853),
               c(66, 0.112, 14.8, 0.112, 0.804))
  # the smallest values are deleted and the remaining excesses over 27.50
  # are kept; the printed data are rounded to one decimal, which explains
  # small differences (k up to 0.014, scale up to 0.22); the first row
  # agrees with the table's W2 = .2309, not with the text's .2389
  s <- sort(y)
  for (i in 1:7) {
    r <- choulakian_stephens2001(s[i:72])
    expect_equal(r$n, pub[i, 1])
    expect_lt(abs(r$k - pub[i, 2]), 0.015)
    expect_lt(abs(r$scale - pub[i, 3]), 0.25)
    expect_lt(abs(r$W2 - pub[i, 4]), 0.005)
    expect_lt(abs(r$A2 - pub[i, 5]), 0.03)
  }
})

test_that("the statistics of Section 2 in closed form, and the p-values", {
  r <- choulakian_stephens2001(y)
  fit <- grimshaw1993(y)
  z <- sort(1 - (1 - fit$k * y / fit$scale)^(1 / fit$k))
  n <- 72
  i <- 1:n
  expect_equal(r$W2, sum((z - (2 * i - 1) / (2 * n))^2) + 1 / (12 * n))
  expect_equal(r$A2, -n - mean((2 * i - 1) * (log(z) + log(1 - rev(z)))))
  # at a tabulated shape, the p-value of a tabulated point is its level
  expect_equal(cs_pvalue(0.974, 0, "A2"), 0.05)
  expect_equal(cs_pvalue(0.153, 0, "W2"), 0.05)
  expect_equal(cs_pvalue(10, 0, "A2"), 0.001)
  expect_equal(cs_pvalue(0.1, 0, "A2"), 0.5)
  # interpolation in k and clamping at 0.5
  expect_equal(cs_points(0.05, "A2"), (cs_points(0, "A2") +
                                         cs_points(0.1, "A2")) / 2)
  expect_equal(cs_points(0.9, "A2"), cs_points(0.5, "A2"))
})

test_that("no estimate, and the exponential case k = 0", {
  expect_error(choulakian_stephens2001(c(1, 2, 3)), "does not exist")
  local_mocked_bindings(grimshaw1993 = function(y) {
    list(exists = TRUE, k = 0, scale = 1)
  })
  r <- choulakian_stephens2001(c(0.1, 0.5, 1, 2, 3))
  z <- 1 - exp(-c(0.1, 0.5, 1, 2, 3))
  i <- 1:5
  expect_equal(r$W2, sum((z - (2 * i - 1) / 10)^2) + 1 / 60)
})

test_that("Table 4: the 5% points of A2 at the estimated shapes", {
  # thresholds 140, 130, 80, 70 of Davison & Smith (1990): k = 0.24, 0.14,
  # -0.34, -0.32 and printed points 1.10, 1.04, 0.865, 0.872
  pts <- vapply(c(0.24, 0.14, -0.34, -0.32), function(k) {
    cs_points(k, "A2")[4]
  }, numeric(1))
  expect_equal(round(pts, 2), c(1.10, 1.04, 0.87, 0.87))
})

test_that("the deletion procedure stops after six deletions", {
  s <- sort(y)
  p <- vapply(0:6, function(d) {
    choulakian_stephens2001(s[(d + 1):72])$p_A2
  }, numeric(1))
  expect_true(all(p[1:6] <= 0.1 + 1e-9))
  expect_gt(p[7], 0.1)
})

test_that("size under a generalized Pareto null", {
  skip_on_cran()
  set.seed(1)
  p <- replicate(1000, {
    x <- (stats::runif(200)^(-0.2) - 1) / 0.2
    choulakian_stephens2001(x)$p_A2
  })
  expect_lt(abs(mean(p < 0.05) - 0.05), 0.02)
})

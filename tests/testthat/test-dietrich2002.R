p8 <- c(.10, .30, .50, .70, .90, .95, .975, .99)
q_at <- function(q, p) {
  stats::approx(evtests:::limit_tables$probs, q, xout = p)$y
}

test_that("the limit quantiles of E reproduce Table 1 of Husler & Li", {
  tab <- evtests:::limit_tables$E
  pub <- list(`0` = c(0.028, 0.042, 0.057, 0.078, 0.122, 0.150, 0.181, 0.222),
              `-0.1` = c(0.027, 0.041, 0.054, 0.074, 0.116, 0.144, 0.174,
                         0.213),
              `-0.3` = c(0.027, 0.040, 0.054, 0.073, 0.113, 0.140, 0.168,
                         0.206),
              `-0.5` = c(0.027, 0.040, 0.054, 0.073, 0.115, 0.141, 0.169,
                         0.208),
              `-0.7` = c(0.028, 0.041, 0.055, 0.074, 0.118, 0.147, 0.176,
                         0.218))
  for (g in names(pub)) {
    q <- q_at(tab$q[, which(abs(tab$gamma - as.numeric(g)) < 1e-9)], p8)
    expect_lt(max(abs(q - pub[[g]]) / pub[[g]]), 0.04)
  }
})

test_that("the limit quantiles of PE reproduce Table 2 of Husler & Li", {
  tab <- evtests:::limit_tables$PE
  pub <- list(`0.5` = c(0.348, 0.493, 0.639, 0.836, 1.250, 1.515, 1.792,
                        2.145),
              `1` = c(0.107, 0.159, 0.211, 0.284, 0.438, 0.534, 0.663, 0.793),
              `2` = c(0.034, 0.053, 0.074, 0.105, 0.176, 0.222, 0.272, 0.341))
  for (e in names(pub)) {
    q <- q_at(tab$q[, match(as.numeric(e), tab$eta)], p8)
    expect_lt(max(abs(q - pub[[e]]) / pub[[e]]), 0.04)
  }
  # the row "> 0" of Table 1 of Dietrich et al. (2002): T = PE with eta = 2
  q <- q_at(tab$q[, 3], c(.90, .95, .975, .99))
  expect_lt(max(abs(q / c(0.173, 0.222, 0.270, 0.330) - 1)), 0.04)
})

test_that("Table 1 of Dietrich et al. for gamma < 0 uses the misprinted sign", {
  # their quantiles (p = .90, .975, .99) are far above those of E with the
  # sign corrected by Husler & Li
  tab <- evtests:::limit_tables$E
  q <- q_at(tab$q[, which(abs(tab$gamma + 0.3) < 1e-9)], c(.90, .975, .99))
  expect_true(all(c(0.396, 0.690, 0.911) / q > 2))
})

test_that("E and PE in closed form", {
  set.seed(1)
  x <- exp(stats::rexp(500))
  z <- sort(x, decreasing = TRUE)
  k <- 30
  d <- log(z[1:k]) - log(z[k + 1])
  m1 <- mean(d)
  m2 <- mean(d^2)
  gm <- 1 - 0.5 / (1 - m1^2 / m2)
  fe <- function(t) {
    j <- pmin(floor(k * t), k - 1)
    (d[j + 1] / m1 - (t^(-gm) - 1) / gm * (1 - gm))^2 * t^2
  }
  fp <- function(t) {
    j <- pmin(floor(k * t), k - 1)
    (d[j + 1] / m1 + log(t))^2 * t^0.5
  }
  piecewise <- function(f) {
    k * sum(vapply(0:(k - 1), function(j) {
      stats::integrate(f, j / k, (j + 1) / k, rel.tol = 1e-12)$value
    }, numeric(1)))
  }
  e <- dietrich2002(x, k = k)
  expect_equal(e$statistic, piecewise(fe), tolerance = 1e-10)
  expect_equal(e$gamma, m1 + gm)
  pe <- dietrich2002(x, k = k, statistic = "PE")
  expect_equal(pe$statistic, piecewise(fp), tolerance = 1e-10)
  # scale invariance (logarithms make them not location invariant)
  expect_equal(dietrich2002(7 * x, k = k)$statistic, e$statistic)
})

test_that("p-values from the tables: the quantile at gamma and the edges", {
  tab <- evtests:::limit_tables$E
  # table_at interpolates linearly between grid columns and holds the edges
  q1 <- tab$q[, which(abs(tab$gamma + 0.1) < 1e-9)]
  q2 <- tab$q[, which(abs(tab$gamma + 0.15) < 1e-9)]
  expect_equal(table_at(tab, -0.125), (q1 + q2) / 2)
  expect_equal(table_at(tab, 0.4), tab$q[, which(tab$gamma == 0)])
  expect_equal(table_at(tab, -2), tab$q[, which.min(tab$gamma)])
  q <- table_at(tab, 0)
  expect_equal(table_pvalue(q_at(q, 0.95), q), 0.05, tolerance = 1e-6)
  expect_equal(table_pvalue(1e6, q), 5e-4)
  expect_equal(table_pvalue(0, q), 1)
})

test_that("size for a Pareto tail and power for a non-EV distribution", {
  skip_on_cran()
  # Dietrich et al. (2002), Figures 3 and 4: Pareto(1) satisfies the
  # condition, exp(floor of an exponential) does not
  set.seed(2)
  p_par <- replicate(400, {
    dietrich2002(top_gpd(10000, 100, 1) + 1, k = 100,
                 statistic = "PE", eta = 2)$p_value
  })
  expect_lt(abs(mean(p_par < 0.05) - 0.05), 0.035)
  p_bad <- replicate(200, {
    x <- exp(floor(-log(stats::runif(10000))))
    suppressWarnings(dietrich2002(x, k = 100, statistic = "PE",
                                  eta = 2)$p_value)
  })
  expect_gt(mean(p_bad < 0.05), 0.8)
})

test_that("validation", {
  expect_error(dietrich2002(1:20, eta = 1), "must be 2")
  expect_error(dietrich2002(1:20, statistic = "PE", eta = 3), "0.5, 1, 2")
  expect_error(dietrich2002(c(5, 4, 3, -1), k = 3), "positive")
  expect_error(dietrich2002(c(1, 1, 1, 0.5), k = 2), "not all be equal")
  r <- dietrich2002(exp(1:20))
  expect_s3_class(r, "evt_test")
  expect_output(print(r), "extreme value condition")
})

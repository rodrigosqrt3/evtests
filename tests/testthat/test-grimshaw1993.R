y <- nylon_fibres$excess

test_that("the summaries of Table 1 (the printed median is a misprint)", {
  expect_equal(round(mean(y), 5), .25267)
  expect_equal(round(stats::sd(y), 5), .24360)
  expect_equal(stats::IQR(y, type = 4), .29)
  expect_equal(stats::median(y), .184)   # printed .16200
})

test_that("Section 4 of Grimshaw (1993) is reproduced", {
  fit <- grimshaw1993(y)
  expect_equal(round(fit$theta_L, 6), -3994.490358)
  expect_equal(round(fit$h2_0, 6), -.008455)
  # no zero on (theta_L, 0); two on (0, theta_U): .415834 (the maximum) and a
  # local minimum close to 1 / max(y), printed as 1.141551, which is not a
  # zero of h (the paper's iteration stopped near the boundary); the zero is
  # at 1.1415244
  z <- fit$zeros
  expect_equal(nrow(z), 2)
  expect_true(all(z$theta > 0))
  # .4158348 here; the paper truncates to .415834
  expect_equal(trunc(z$theta[1] * 1e6) / 1e6, .415834)
  expect_equal(round(z$theta[2], 7), 1.1415244)
  expect_equal(z$maximum, c(TRUE, FALSE))
  expect_equal(round(z$k[1], 6), .117698)
  expect_equal(round(z$alpha[1], 6), .283040)
  expect_equal(round(z$loglik[1], 6), 5.697968)
  # the boundary value -n log(max(y)) and the decision of step 8
  expect_equal(round(-15 * log(max(y)), 6), 1.985838)
  expect_true(fit$exists)
  expect_false(fit$boundary)
  expect_equal(round(fit$k, 6), .117698)
  expect_equal(round(fit$scale, 6), .283040)
  expect_equal(fit$gamma, -fit$k)
  # Q(1) = alpha / k = 2.404804 printed; 2.404801 at full precision
  expect_equal(round(fit$endpoint, 5), 2.40480)
  # "approximate 95% upper confidence bound" 9.099524: one-sided, z = 1.645
  ub <- fit$endpoint + stats::qnorm(0.95) * fit$endpoint_se
  expect_equal(round(ub, 3), 9.099)
})

test_that("the estimate solves the likelihood equations", {
  fit <- grimshaw1993(y)
  k <- fit$k
  a <- fit$scale
  n <- length(y)
  # the gradient of the log-likelihood (appendix) vanishes
  w <- 1 / (1 - k * y / a)
  dk <- (n / k) * (1 / k - 1) - sum(log(1 - k * y / a)) / k^2 -
    (1 / k) * (1 / k - 1) * sum(w)
  da <- -n / (k * a) + (1 / a) * (1 / k - 1) * sum(w)
  expect_lt(abs(dk), 1e-8)
  expect_lt(abs(da), 1e-8)
  # the same equations in the notation of Drees, Ferreira & de Haan (2004)
  g <- fit$gamma
  expect_equal(mean(log(1 + g / a * y)), g)
  expect_equal(mean(1 / (1 + g / a * y)), 1 / (g + 1))
  # and it is the maximum of the log-likelihood
  ll <- function(p) -n * log(p[2]) + (1 / p[1] - 1) * sum(log(1 - p[1] * y / p[2]))
  opt <- stats::optim(c(0.2, 0.3), function(p) {
    if (p[2] <= 0 || any(1 - p[1] * y / p[2] <= 0)) return(1e10)
    -ll(p)
  }, control = list(reltol = 1e-14))
  expect_equal(fit$loglik, -opt$value, tolerance = 1e-8)
  expect_equal(fit$loglik, ll(c(k, a)))
})

test_that("asymptotic covariance of Smith (1984) as quoted", {
  fit <- grimshaw1993(y)
  k <- fit$k
  a <- fit$scale
  expect_equal(unname(fit$vcov),
               matrix(c((1 - k)^2, a * (1 - k), a * (1 - k), 2 * a^2 * (1 - k)),
                      2) / 15)
  expect_equal(fit$endpoint_se^2, a^2 * (1 - k)^2 * (1 - 2 * k) / (15 * k^4))
})

test_that("heavy, exponential and short tails", {
  set.seed(1)
  heavy <- grimshaw1993(stats::runif(2000)^(-0.5) - 1)  # GPD with gamma = 1/2
  expect_equal(heavy$gamma, 0.5, tolerance = 0.1)
  expect_lt(heavy$theta, 0)
  expect_true(is.na(heavy$endpoint))
  ex <- grimshaw1993(stats::rexp(2000, rate = 1 / 3))
  expect_lt(abs(ex$gamma), 0.08)
  expect_equal(ex$scale, 3, tolerance = 0.1)
  short <- grimshaw1993(1 - stats::runif(2000)^(1 / 3)) # gamma = -1/3
  expect_equal(short$gamma, -1 / 3, tolerance = 0.1)
  expect_gt(short$endpoint, max(1 - stats::runif(1)))
  expect_true(is.finite(short$endpoint_se))
})

test_that("the boundary point and a sample without an estimate", {
  # a local maximum (k = .511) whose likelihood is below the boundary value
  # -n log(max(u)), so that step 8 chooses (k, alpha) = (1, max(u))
  u <- c(0.659, 0.199, 0.05, 0.347, 0.249, 0.655, 0.273, 0.034, 0.343, 0.155,
         0.211, 0.163, 0.002)
  b <- grimshaw1993(u)
  expect_true(b$zeros$maximum[1])
  expect_lt(b$zeros$loglik[1], -13 * log(max(u)))
  expect_true(b$boundary)
  expect_equal(b$k, 1)
  expect_equal(b$scale, max(u))
  expect_equal(b$loglik, -13 * log(max(u)))
  expect_true(all(is.na(b$vcov)))
  expect_output(print(b), "boundary")
  # without a local maximum there is no estimate (step 7), even though the
  # boundary point exists: uniform excesses
  set.seed(2)
  expect_false(grimshaw1993(stats::runif(50))$exists)
  # equally spaced excesses have no stationary maximum with k <= 1
  none <- grimshaw1993(c(1, 2, 3))
  expect_false(none$exists)
  expect_true(is.na(none$gamma))
  expect_output(print(none), "no maximum")
  expect_output(print(grimshaw1993(y)), "Grimshaw")
})

test_that("input validation", {
  expect_error(grimshaw1993("a"), "numeric")
  expect_error(grimshaw1993(c(1, NA)), "numeric")
  expect_error(grimshaw1993(c(-1, 2)), "non-negative")
  expect_error(grimshaw1993(c(0, 0)), "non-negative")
  expect_error(grimshaw1993(1), "at least 2")
})

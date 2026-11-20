test_that("the estimate is the end-point of the generalized Pareto MLE", {
  set.seed(1)
  x <- 1 - stats::runif(500)^(1 / 3)
  k <- 60
  s <- smith1985(x, k = k)
  z <- sort(x, decreasing = TRUE)
  g <- grimshaw1993(z[1:k] - z[k + 1])
  expect_gt(g$k, 0)
  expect_equal(s$endpoint, z[k + 1] + g$scale / g$k, tolerance = 1e-6)
  expect_equal(s$gamma, g$gamma, tolerance = 1e-5)
  expect_equal(s$shape, -1 / g$gamma, tolerance = 1e-5)
  expect_true(s$regular)
})

test_that("the profile log-likelihood and the interval", {
  set.seed(2)
  x <- 1 - stats::runif(800)^(1 / 4)
  k <- 80
  s <- smith1985(x, k = k, alpha = 0.1)
  z <- sort(x, decreasing = TRUE)
  top <- z[1:k]
  u <- z[k + 1]
  # the profile is the Beta(1, a) log-likelihood maximised over a
  th <- s$endpoint
  r <- (th - top) / (th - u)
  a <- -k / sum(log(r))
  ll <- sum(stats::dbeta((top - u) / (th - u), 1, a, log = TRUE)) -
    k * log(th - u)
  expect_equal(s$loglik, ll, tolerance = 1e-8)
  # the limits are where the deviance equals the chi-squared quantile
  for (lim in c(s$lower, s$upper)) {
    expect_equal(2 * (s$loglik - smith_profile(lim - top[1], top, u)),
                 stats::qchisq(0.9, 1), tolerance = 1e-6)
  }
  expect_lt(s$lower, s$endpoint)
  expect_gt(s$upper, s$endpoint)
})

test_that("coverage for shape 4 and a warning in the non-regular case", {
  skip_on_cran()
  set.seed(3)
  cov <- replicate(300, {
    x <- 1 - stats::runif(1000)^(1 / 4)
    r <- suppressWarnings(smith1985(x, k = 100))
    isTRUE(r$lower <= 1 && 1 <= r$upper)
  })
  expect_gt(mean(cov), 0.85)
  # shape 1 (uniform): the estimate usually does not exist, and when it does
  # the fitted shape is at most about 2 and a warning is given
  set.seed(4)
  out <- vapply(1:20, function(i) {
    w <- NULL
    r <- withCallingHandlers(smith1985(stats::runif(300), k = 60),
                             warning = function(c) {
                               w <<- conditionMessage(c)
                               invokeRestart("muffleWarning")
                             })
    is.na(r$endpoint) || (!is.null(w) && grepl("non-regular", w)) ||
      r$shape > 2
  }, logical(1))
  expect_true(all(out))
  expect_gt(mean(vapply(1:20, function(i) {
    is.na(suppressWarnings(smith1985(stats::runif(300), k = 60))$endpoint)
  }, logical(1))), 0.5)
})

test_that("no interior maximum, the non-regular warning, a lower limit of 0", {
  na <- smith1985(c(1.1, 0.8, 1.7, 8.4, 5.4, 5.1, 7.6, 2.5), k = 6)
  expect_true(is.na(na$endpoint) && is.na(na$lower))
  expect_warning(smith1985(c(0.9, 1.2, 5.1, 4.9, 2.7, 3.5, 0.6, 10), k = 6),
                 "non-regular")
  lo <- suppressWarnings(smith1985(c(0.4, 3.3, 6.9, 6.3, 3.4, 9.3, 3.2, 3.8),
                                   k = 6))
  expect_equal(lo$lower, 9.3)
  expect_equal(suppressWarnings(smith1985(stats::rexp(30)))$k, 10L)
})

test_that("an interior maximum below the limit at infinity is not kept", {
  set.seed(5)
  x <- sort(1 - stats::runif(200)^(1 / 3), decreasing = TRUE)
  top <- x[1:30]
  u <- x[31]
  expect_false(is.na(smith_fit(top, u, top[1] - u)$s))
  local_mocked_bindings(smith_profile_inf = function(x, u) Inf)
  expect_true(is.na(smith_fit(top, u, top[1] - u)$s))
})

test_that("validation", {
  expect_error(smith1985(1:4), "at least 5")
  expect_error(smith1985(1:20, k = 3), "between 4")
  expect_error(smith1985(1:20, k = c(5, 6)), "single")
  expect_error(smith1985(c(rep(2, 6), 1), k = 5), "not all be equal")
})

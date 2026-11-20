test_that("Y, the Hill path, T (2.3) and Q (2.6) in closed form", {
  set.seed(1)
  x <- (-log(stats::runif(300)))^(-1)
  res <- guillou_hall2001(x)
  z <- sort(x, decreasing = TRUE)
  k <- 40
  y <- (1:k) * (log(z[1:k]) - log(z[2:(k + 1)]))
  expect_equal(res$path$gamma[k], hill1975(x, k = k)$gamma)
  w <- sign(k - 2 * (1:k) + 1) * abs(k - 2 * (1:k) + 1)
  expect_equal(res$path$T[k], sqrt(3 / k^3) * sum(w * y) / mean(y))
  l <- k %/% 2
  expect_equal(res$path$Q[k], sqrt(mean(res$path$T[(k - l):(k + l)]^2)))
  # the rule "first"
  q <- res$path$Q
  expect_equal(res$k, which(!is.na(q) & q >= 1.25 & seq_along(q) >= 2)[1])
  expect_equal(res$gamma, res$path$gamma[res$k])
})

test_that("T is approximately standard normal for an exact Pareto tail", {
  set.seed(2)
  t50 <- replicate(3000, guillou_hall2001(stats::runif(200)^(-0.5),
                                          statistic = "T")$path$T[50])
  expect_lt(abs(mean(t50)), 0.06)
  expect_lt(abs(stats::sd(t50) - 1), 0.05)
})

test_that("Table 1: mean selected k for Frechet and Burr, n = 500", {
  skip_on_cran()
  set.seed(3)
  mk <- function(gen) mean(replicate(800, guillou_hall2001(gen())$k))
  # printed (crit = 1.25): Frechet(1) 90, Burr(1, 1, 1) 69
  expect_lt(abs(mk(function() (-log(stats::runif(500)))^(-1)) - 90), 10)
  expect_lt(abs(mk(function() 1 / stats::runif(500) - 1) - 69), 8)
})

test_that("the rule 'all', T, and validation", {
  set.seed(4)
  x <- (-log(stats::runif(500)))^(-1)
  a <- guillou_hall2001(x, rule = "all")
  q <- a$path$Q
  ok <- which(!is.na(q))
  expect_true(all(q[ok][ok >= a$k] > 1.25))
  expect_true(q[a$k - 1] <= 1.25)
  tt <- guillou_hall2001(x, statistic = "T", crit = 1.5)
  expect_true(abs(tt$path$T[tt$k]) >= 1.5)
  expect_true(is.na(guillou_hall2001(x, crit = 100)$k))
  expect_error(guillou_hall2001(x, crit = -1), "positive")
  expect_error(guillou_hall2001(c(-1, -2, 1, 2)), "4 positive")
})

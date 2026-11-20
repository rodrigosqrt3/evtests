test_that("Corollary 7 equals (3.10) with the rho of (3.9) and Corollary 3", {
  for (case in list(c(1000, 40, 9), c(1500, 80, 30))) {
    n1 <- case[1]
    k1 <- case[2]
    k2 <- case[3]
    cor7 <- (k1^2 / k2) *
      ((log(k1))^2 / (2 * log(n1) - log(k1))^2)^((log(n1) - log(k1)) /
                                                    log(n1))
    rho <- log(k1) / (-2 * log(n1) + 2 * log(k1))
    expect_equal(cor7, (k1^2 / k2) / (1 - 1 / rho)^(2 / (1 - 2 * rho)))
  }
})

test_that("the bootstrap criterion is M - 2 gamma^2 averaged over resamples", {
  set.seed(1)
  z <- sort(stats::runif(50)^(-0.5), decreasing = TRUE)
  set.seed(2)
  q <- dan_q(z, 30, 5)
  set.seed(2)
  vals <- replicate(5, {
    s <- sort(sample(z, 30, replace = TRUE), decreasing = TRUE)
    vapply(1:29, function(k) {
      e <- log(s[1:k]) - log(s[k + 1])
      (mean(e^2) - 2 * mean(e)^2)^2
    }, numeric(1))
  })
  crit <- rowMeans(vals)
  expect_equal(q$k, which.min(crit))
  expect_equal(q$q, min(crit))
})

test_that("the selection over n1 and the returned estimate", {
  set.seed(3)
  x <- (-log(stats::runif(1000)))^(-0.5)
  res <- danielsson2001(x, n1 = c(400, 600), B = 50)
  g <- res$grid
  expect_equal(g$n2, round(c(400, 600)^2 / 1000))
  expect_equal(res$n1, g$n1[which.min(g$R)])
  expect_equal(res$gamma, hill1975(x, k = res$k)$gamma)
  expect_output(print(res), "selected k")
})

test_that("Table I, Extreme(4), n = 2000: mean gamma .259", {
  skip_on_cran()
  set.seed(4)
  g <- replicate(30, danielsson2001((-log(stats::runif(2000)))^(-0.25),
                                    B = 150)$gamma)
  # printed mean .259, standard error .024 over 250 simulations
  expect_lt(abs(mean(g) - 0.259), 3 * 0.024 / sqrt(30) + 0.005)
})

test_that("validation", {
  expect_error(danielsson2001(1:5), "at least 10")
  expect_error(danielsson2001(stats::runif(50), n1 = 60), "between 2")
  expect_error(danielsson2001(stats::runif(50), B = 0), "positive whole")
})

frechet <- function(n, beta) (-log(stats::runif(n)))^(-beta)

test_that("the stopping time, the initial estimate and k_opt in closed form", {
  set.seed(1)
  x <- frechet(1000, 0.5)
  res <- drees_kaufmann1998(x)
  z <- sort(x, decreasing = TRUE)
  b <- function(i) mean(log(z[1:(i - 1)])) - log(z[i])
  g0 <- b(floor(2 * sqrt(1000)))
  expect_equal(res$gamma_initial, g0)
  kbar <- function(r) {
    repeat {
      for (k in 2:1000) {
        if (max(sqrt(2:k) * abs(vapply(2:k, b, 0) - b(k))) > r) {
          return(c(k, r))
        }
      }
      r <- 0.9 * r
    }
  }
  r <- 2.5 * g0 * 1000^0.25
  k1 <- kbar(r)
  k2 <- kbar(r^0.7)
  expect_equal(unname(res$kbar), c(k1[1], k2[1]))
  kopt <- floor(3^(-1) * (2 * g0^2)^(1 / 3) * (k2[1] / k1[1]^0.7)^(1 / 0.3))
  expect_equal(res$k, kopt - 1)
  expect_equal(res$gamma, b(kopt))
  expect_equal(res$gamma, hill1975(x, k = res$k)$gamma)
  expect_output(print(res), "selected k")
})

test_that("Tables 1-4 of Drees & Kaufmann: Frechet, n = 500 and 1000", {
  skip_on_cran()
  sim <- function(beta, n, reps) {
    out <- t(replicate(reps, {
      r <- drees_kaufmann1998(frechet(n, beta))
      c(r$k, r$gamma)
    }))
    ok <- !is.na(out[, 1])
    c(fail = mean(!ok), rmse = sqrt(mean((out[ok, 2] - beta)^2)),
      kmed = stats::median(out[ok, 1] + 1))
  }
  set.seed(2)
  # printed: F0.5, n = 1000: 0.0% failures, RMSE .048, median k / k_opt .86
  # (k_opt = 2 n^(2/3) = 200, counting the threshold)
  a <- sim(0.5, 1000, 600)
  expect_lt(a[["fail"]], 0.01)
  expect_lt(abs(a[["rmse"]] - 0.048), 0.005)
  expect_lt(abs(a[["kmed"]] / 200 - 0.86), 0.12)
  # F0.2, n = 500: RMSE .027, median k / k_opt .5 (k_opt = 126)
  b <- sim(0.2, 500, 600)
  expect_lt(abs(b[["rmse"]] - 0.027), 0.004)
  expect_lt(abs(b[["kmed"]] / 126 - 0.5), 0.1)
})

test_that("estimated rho, other constants and validation", {
  set.seed(3)
  x <- frechet(2000, 1)
  e <- drees_kaufmann1998(x, rho = "estimate")
  expect_true(e$rho > 0)
  expect_true(!is.na(e$k))
  expect_equal(dk_rho1(c(NA, rep(1, 10)), 4, 0.6), 1)
  # a flat Hill path gives 0 / 0, so rho falls back to 1
  expect_equal(dk_rho1(c(NA, rep(1, 10)), 10, 0.6), 1)
  # a threshold too large for the path is reduced by 0.9 until kbar exists
  big <- drees_kaufmann1998(x, r_factor = 1e3)
  expect_lt(big$r[["r"]], 1e3 * big$gamma_initial * 2000^0.25)
  # tied top observations give a zero initial estimate
  expect_error(drees_kaufmann1998(c(rep(5, 200), 1:3)), "not positive")
  expect_error(drees_kaufmann1998(x, rho = -1), "positive")
  expect_error(drees_kaufmann1998(x, xi = 1), "xi")
  expect_error(drees_kaufmann1998(c(-1, -2, 3, 4)), "4 positive")
})

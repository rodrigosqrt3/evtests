test_that("ForwardStop (equation 2) and the selected threshold", {
  p <- c(0.001, 0.01, 0.2, 0.5, 0.9)
  fs <- forward_stop(p)
  expect_equal(fs, -cumsum(log(1 - p)) / (1:5))
  # by hand: kF is the largest j with fs <= alpha
  expect_equal(max(which(fs <= 0.05)), 2L)
  expect_equal(max(which(fs <= 0.2)), 3L)
})

test_that("bader2018 applies ForwardStop to the Choulakian-Stephens p-values", {
  set.seed(1)
  x <- c(stats::rexp(400, 3), 1 + stats::rexp(300))
  res <- bader2018(x, k = c(600, 400, 300, 200, 100, 50))
  z <- sort(x, decreasing = TRUE)
  p <- vapply(res$path$k, function(kk) {
    choulakian_stephens2001(z[1:kk] - z[kk + 1])$p_A2
  }, 0)
  expect_equal(res$path$p_value, p)
  expect_equal(res$path$k, c(600, 400, 300, 200, 100, 50))
  kf <- max(c(0, which(forward_stop(p) <= 0.05)))
  expect_equal(res$kF, kf)
  expect_equal(res$k, res$path$k[if (kf == 6) 6 else kf + 1])
  expect_equal(res$threshold, z[res$k + 1])
  expect_output(print(res), "selected k")
})

test_that("no rejection selects the lowest threshold; all, the highest", {
  set.seed(2)
  g <- (stats::runif(800)^(-0.1) - 1) / 0.1
  r <- bader2018(g, k = c(700, 500, 300))
  if (r$kF == 0) expect_equal(r$k, 700)
  set.seed(3)
  w <- stats::rweibull(2000, 0.5)
  rw <- bader2018(w, k = c(1500, 1200, 1000))
  if (rw$kF == 3) expect_equal(rw$k, 1000)
  expect_error(bader2018(1:20, k = 5), "between 10")
})

test_that("equations 3 and 7, and the Gumbel critical values", {
  set.seed(1)
  x <- stats::rlnorm(300)
  z <- sort(x, decreasing = TRUE)
  k <- 40
  e <- z[1:k] - z[k + 1]
  res <- neves2006(x, k = k, alternative = "positive")
  expect_equal(res$statistic, e[1] / mean(e) - log(k))
  expect_equal(res$p_value, 1 - exp(-exp(-res$statistic)))
  # the two-sided 5% rejection region of p. 1286: below g(.025) = -1.3053 or
  # above g(.975) = 3.6762, g(e) = -log(-log(e))
  expect_equal(round(-log(-log(c(.025, .975))), 4), c(-1.3053, 3.6762))
  expect_equal(neves2006(5 + 2 * x, k = k)$statistic, res$statistic)
})

test_that("the exact null for an exponential tail (Fisher's g)", {
  # P(max / sum > x) by simulation of k exponentials
  set.seed(2)
  for (k in c(5, 40, 300)) {
    e <- matrix(stats::rexp(2e4 * k), ncol = k)
    r <- apply(e, 1, max) / rowSums(e)
    x <- stats::quantile(r, c(.1, .5, .95), names = FALSE)
    p <- vapply(x, fisher_g_upper, numeric(1), k = k)
    expect_lt(max(abs(p - c(.9, .5, .05))), 0.012)
  }
  expect_equal(fisher_g_upper(1, 10), 0)
  expect_equal(fisher_g_upper(0.05, 10), 1)
  # far in the lower tail for large k the terms overflow and NA is returned
  expect_true(is.na(fisher_g_upper(1.5 / 5000, 5000)))
})

test_that("size under an exponential tail and consistency (Theorems 6, 9)", {
  skip_on_cran()
  set.seed(3)
  p <- replicate(2000, neves2006(top_gpd(1000, 50, 0), k = 50)$p_exact)
  expect_lt(abs(mean(p < 0.05) - 0.05), 0.015)
  expect_lt(neves2006(stats::runif(2000), k = 200,
                      alternative = "negative")$p_value, 0.001)
  expect_lt(neves2006(stats::runif(2000)^(-1), k = 200,
                      alternative = "positive")$p_value, 0.001)
})

test_that("structure and validation", {
  set.seed(4)
  res <- neves2006(stats::rexp(30))
  expect_s3_class(res, "evt_test")
  expect_equal(res$k, 2:29)
  expect_true(all(c("p_value", "p_exact") %in% names(res)))
  expect_error(neves2006(c(1, 1, 1, 0), k = 2), "not all be equal")
})

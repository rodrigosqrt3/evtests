city <- us_cities_1940$population

test_that("Table 1 of Hill (1975): iV, alpha0 and alpha1 to the printed digits", {
  h <- hill1975(city, k = 1:29)
  iv_pub <- c(.79, 1.13, .52, .30, 2.69, .13, .36, .45, 1.24, .13, .47, .94,
              .25, 2.12, .09, 1.22, 1.00, 1.35, .58, .10, .95, 2.73, .21, .22,
              1.04, .09, .24, .01, .10)
  a0_pub <- c(1.27, 1.04, 1.23, 1.46, .92, 1.08, 1.18, 1.25, 1.18, 1.29, 1.34,
              1.31, 1.38, 1.21, 1.29, 1.25, 1.23, 1.19, 1.21, 1.26, 1.25, 1.13,
              1.16, 1.20, 1.19, 1.23, 1.26, 1.31, 1.35)
  a1_pub <- c(2.54, 1.57, 1.64, 1.82, 1.10, 1.26, 1.35, 1.41, 1.31, 1.42, 1.46,
              1.42, 1.49, 1.30, 1.38, 1.32, 1.30, 1.25, 1.27, 1.32, 1.31, 1.18,
              1.21, 1.25, 1.24, 1.28, 1.31, 1.36, 1.40)
  expect_equal(round(attr(h, "iV"), 2), iv_pub)
  expect_equal(round(h$alpha0, 2), a0_pub)
  expect_equal(round(h$alpha1, 2), a1_pub)
  expect_equal(h$gamma, 1 / h$alpha0)
})

test_that("Table 1 of Hill (1975): the discrepancy statistic H", {
  h <- hill1975(city, k = 2:29)
  h_pub <- c(.06, .28, .81, 3.05, 5.07, 6.47, 7.51, 6.90, 8.99, 9.80, 9.45,
             10.99, 11.16, 13.45, 12.82, 12.52, 12.03, 12.53, 14.50, 14.27,
             16.10, 17.81, 19.52, 19.19, 21.42, 23.09, 25.81, 28.23)
  # 24 of the 28 printed values agree at two decimals; k = 15, 26, 28, 29
  # differ by 0.005 to 0.015 (13.445, 21.425, 25.825, 28.239 here), which the
  # rounded iV of the table do not explain either
  expect_lt(max(abs(h$H - h_pub)), 0.02)
  expect_equal(sum(round(h$H, 2) == h_pub), 24)
  expect_equal(hill1975(city, k = 1)$H, 0)
})

test_that("Table 1 of Hill (1975): k C1, in inhabitants", {
  h <- hill1975(city, k = 1:29)
  kc_pub <- c(8.3e16, 2.1e10, 6.3e10, 9.2e11, 2.2e7, 2.0e8, 7.7e8, 1.8e9,
              4.5e8, 2.0e9, 3.5e9, 2.0e9, 5.2e9, 3.8e8, 1.1e9, 5.3e8, 3.8e8,
              1.9e8, 2.4e8, 5.2e8, 4.2e8, 7.0e7, 1.2e8, 1.9e8, 1.6e8, 2.8e8,
              4.3e8, 7.9e8, 1.4e9)
  kc <- signif(h$C1 * length(city), 2)
  # two rows differ in the second digit (k = 3: 6.24e10, k = 28: 7.97e8), very
  # sensitive to alpha1, which enters as an exponent of about 15
  expect_equal(which(kc != kc_pub), c(3L, 28L))
  expect_equal(kc[-c(3, 28)], kc_pub[-c(3, 28)])
})

test_that("the text results tied to Table 1", {
  h <- hill1975(city, k = 29)
  # "the standard deviation of alpha0(29), as given by (2.6), is .20 alpha"
  expect_equal(round(h$sd_rel, 2), 0.20)
  # chi-square of the 29 iV against an exponential, classes (0, .5], (.5, 1]
  # and above 1: "a value of .3" (the paper does not say which rate it used;
  # both estimates give .3)
  iv <- attr(h, "iV")
  obs <- table(cut(iv, c(0, 0.5, 1, Inf)))
  chisq <- function(rate) {
    p <- diff(stats::pexp(c(0, 0.5, 1, Inf), rate))
    sum((obs - 29 * p)^2 / (29 * p))
  }
  expect_equal(as.vector(obs), c(15, 6, 8))
  expect_equal(round(chisq(h$alpha0), 1), 0.3)
  expect_equal(round(chisq(h$alpha1), 1), 0.3)
})

test_that("alpha0 has the conditional law of equation 2.6 under a Pareto tail", {
  set.seed(1)
  k <- 20
  a0 <- replicate(4000, hill1975(stats::runif(200)^(-1 / 2), k = k)$alpha0)
  # E = k alpha / (k - 1); sd = alpha k / ((k - 1) sqrt(k - 2))
  expect_equal(mean(a0), 2 * k / (k - 1), tolerance = 0.02)
  expect_equal(stats::sd(a0), 2 * hill1975(1:30, k = k)$sd_rel,
               tolerance = 0.05)
})

test_that("structure, the whole path, printing and plotting", {
  h <- hill1975(city)
  expect_s3_class(h, "evt_estimate")
  expect_equal(h$k, 1:29)
  expect_true(all(is.na(h$sd_rel[1:2])))
  expect_output(print(hill1975(city, k = 5)), "Hill \\(1975\\)")
  pdf(NULL)
  on.exit(dev.off())
  expect_invisible(plot(h))
  expect_error(plot(h, "nothing"), "no column")
})

test_that("input validation", {
  expect_error(hill1975("a"), "numeric")
  expect_error(hill1975(c(1, NA, 3)), "missing")
  expect_error(hill1975(1), "at least 2")
  expect_error(hill1975(city, k = 30), "between 1 and 29")
  expect_error(hill1975(city, k = 2.5), "whole")
  expect_error(hill1975(c(5, 4, 0, -1), k = 3), "positive")
  expect_equal(nrow(hill1975(city, k = c(5, 3, 5))), 2)
})

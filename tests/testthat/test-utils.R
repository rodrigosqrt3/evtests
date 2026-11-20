test_that("evt_k validates k and refuses an empty range", {
  expect_equal(evt_k(NULL, 2L, 5L, "m"), 2:5)
  expect_equal(evt_k(c(4, 2, 4), 2L, 5L, "m"), c(2L, 4L))
  expect_error(evt_k(NULL, 3L, 2L, "m"), "admissible k range is empty")
  expect_error(evt_k(1, 2L, 5L, "m"), "between 2 and 5")
  expect_error(evt_k(2.5, 2L, 5L, "m"), "whole")
  # for example the bound of Fraga Alves et al. needs k >= 2, n / 2 >= 2
  expect_error(fraga_alves_neves2014(1:3, bound = TRUE), "empty")
})

test_that("evt_sorted and check_alpha", {
  expect_equal(evt_sorted(c(2, 3, 1), "m"), c(3, 2, 1))
  expect_error(evt_sorted(c(1, Inf), "m"), "infinite")
  expect_error(check_alpha(c(0.1, 0.2)), "single")
})

test_that("sim_null caches, uses a fixed seed and leaves no seed behind", {
  had <- exists(".Random.seed", envir = globalenv(), inherits = FALSE)
  if (had) old <- get(".Random.seed", envir = globalenv())
  on.exit(if (had) assign(".Random.seed", old, envir = globalenv()))
  rm(".Random.seed", envir = globalenv())
  a <- sim_null("unit-test-a", 200, function(ns) stats::runif(ns))
  expect_false(exists(".Random.seed", envir = globalenv(), inherits = FALSE))
  expect_equal(a, sort(a))
  # cached: the generator is not called again
  b <- sim_null("unit-test-a", 200, function(ns) stop("not cached"))
  expect_identical(a, b)
  # the same seed gives the same draws under another key
  set.seed(20261120L)
  expect_equal(sim_null("unit-test-b", 200, function(ns) stats::runif(ns)),
               sort(stats::runif(200)))
})

test_that("evt_pick and the alternatives", {
  expect_equal(evt_pick(0.2, 0.9, "positive"), 0.2)
  expect_equal(evt_pick(0.2, 0.9, "negative"), 0.9)
  expect_equal(evt_pick(0.2, 0.9, "two.sided"), 0.4)
  expect_equal(evt_pick(0.7, 0.8, "two.sided"), 1)
  expect_error(evt_alternative("up"), "should be one")
})

test_that("a fixed nsims is still honoured exactly", {
  res <- suppressMessages(
    adsim(ngens = 20, k = 1, l = 20, admix = 0.5, m1 = 0, m2 = 0, ne = 200,
          p1 = 0.5, p2 = 0.5, nsims = 300, keep_sims = FALSE)
  )
  expect_equal(res$conditions$nsims, 300)
  expect_equal(res$summary$nsims, 300)
  # the precision reached is reported even when it was not targeted
  expect_true(is.finite(res$summary$mc_precision))
})

test_that("auto reaches the precision it was asked for", {
  res <- suppressMessages(
    adsim(ngens = 33, k = 1.1, l = 5, admix = 0.5, m1 = 0, m2 = 0, ne = 200,
          p1 = 0.5, p2 = 0.9, nsims = "auto", precision = 0.02,
          keep_sims = FALSE)
  )
  expect_true(res$summary$converged)
  expect_lte(res$summary$mc_precision, 0.02)
  expect_gte(res$conditions$nsims, 500)
})

test_that("a tighter precision demands more runs", {
  params <- list(ngens = 33, k = 1.1, l = 5, admix = 0.5, m1 = 0, m2 = 0,
                 ne = 200, p1 = 0.5, p2 = 0.9, nsims = "auto",
                 keep_sims = FALSE)
  loose <- suppressMessages(do.call(adsim, c(params, list(precision = 0.05))))
  tight <- suppressMessages(do.call(adsim, c(params, list(precision = 0.01))))
  expect_gt(tight$conditions$nsims, loose$conditions$nsims)
  expect_lte(tight$summary$mc_precision, 0.01)
})

test_that("each scenario gets the number of runs it needs", {
  # A distribution concentrated by migration in a large population is resolved
  # by far fewer runs than one spread across the whole interval.
  easy <- suppressMessages(
    adsim(ngens = 20, k = 1, l = 20, admix = 0.5, m1 = 0.05, m2 = 0.05,
          ne = 5000, p1 = 0.5, p2 = 0.5, nsims = "auto", keep_sims = FALSE)
  )
  hard <- suppressMessages(
    adsim(ngens = 33, k = 1.1, l = 5, admix = 0.5, m1 = 0, m2 = 0, ne = 200,
          p1 = 0.5, p2 = 0.9, nsims = "auto", keep_sims = FALSE)
  )
  expect_lt(easy$conditions$nsims, hard$conditions$nsims)
})

test_that("max_nsims stops the search and is reported as not converged", {
  res <- suppressMessages(
    adsim(ngens = 33, k = 1.1, l = 5, admix = 0.5, m1 = 0, m2 = 0, ne = 200,
          p1 = 0.5, p2 = 0.9, nsims = "auto", precision = 1e-4,
          min_nsims = 200, max_nsims = 600, keep_sims = FALSE)
  )
  expect_false(res$summary$converged)
  expect_lte(res$conditions$nsims, 600 + 200)
  expect_gt(res$summary$mc_precision, 1e-4)
})

test_that("a quantile stranded between two atoms is excluded and counted", {
  # Half the runs lost and half fixed puts the median exactly on the jump, so
  # its interval spans 0 to 1 whatever the number of runs. Levels buried
  # inside an atom are determined and must still count.
  split <- c(rep(0, 500), rep(1, 500))
  got <- run_precision(split, precision = 0.01)
  expect_equal(got$at_atom, 1L)
  expect_equal(got$half_width_quantiles, 0)
  expect_equal(got$half_width, got$half_width_props)

  # every run at the same value: nothing is left undetermined
  single <- run_precision(rep(0.3, 500), precision = 0.01)
  expect_equal(single$at_atom, 0L)
  expect_equal(single$half_width_quantiles, 0)

  # a spread out distribution has no atoms at all
  spread <- run_precision(seq(0, 1, length.out = 2000), precision = 0.01)
  expect_equal(spread$at_atom, 0L)
  expect_gt(spread$half_width_quantiles, 0)
})

test_that("an absorbed scenario converges on its proportions", {
  # Runs taken past 4N are all lost or fixed, so the distribution is made of
  # its two atoms and the proportions are what must converge.
  res <- suppressMessages(
    adsim(ngens = 1500, k = 1, l = 1500, admix = 0.5, m1 = 0, m2 = 0, ne = 50,
          p1 = 0.5, p2 = 0.5, nsims = "auto", precision = 0.02,
          fixed_p0 = 0.5, keep_sims = FALSE)
  )
  expect_equal(res$summary$prop_segregating, 0)
  expect_true(res$summary$converged)
  expect_lte(res$summary$mc_precision_props, 0.02)
})

test_that("precision arguments are validated", {
  base <- list(ngens = 10, k = 1, l = 10, admix = 0.5, m1 = 0, m2 = 0,
               ne = 100, p1 = 0.5, p2 = 0.5)
  expect_error(do.call(adsim, c(base, list(nsims = "auto", precision = 0))))
  expect_error(do.call(adsim, c(base, list(nsims = "auto", precision = 2))))
  expect_error(do.call(adsim, c(base, list(nsims = "auto", min_nsims = 1))))
  expect_error(do.call(adsim, c(base, list(nsims = -5))))
  expect_error(do.call(adsim, c(base, list(nsims = 10.5))))
})

test_that("Wilson and order statistic half widths behave", {
  # a proportion known exactly still carries the usual small-sample width
  expect_lt(wilson_half_width(0, 1000), 0.01)
  expect_gt(wilson_half_width(500, 1000), 0.02)
  # a degenerate sample is resolved immediately
  flat <- run_precision(rep(0.5, 1000), precision = 0.01)
  expect_equal(flat$half_width_quantiles, 0)
})

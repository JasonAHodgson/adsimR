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
  # tail_p is relaxed so that the shape criterion, not the tail depth, is what
  # decides the number of runs here.
  params <- list(ngens = 33, k = 1.1, l = 5, admix = 0.5, m1 = 0, m2 = 0,
                 ne = 200, p1 = 0.5, p2 = 0.9, nsims = "auto",
                 tail_p = 0.2, keep_sims = FALSE)
  loose <- suppressMessages(do.call(adsim, c(params, list(precision = 0.05))))
  tight <- suppressMessages(do.call(adsim, c(params, list(precision = 0.01))))
  expect_gt(tight$conditions$nsims, loose$conditions$nsims)
  expect_lte(tight$summary$mc_precision, 0.01)
})

test_that("each scenario gets the number of runs it needs", {
  # A distribution concentrated by migration in a large population is resolved
  # by far fewer runs than one spread across the whole interval.
  # With the tail depth relaxed, the shape of the distribution decides.
  easy <- suppressMessages(
    adsim(ngens = 20, k = 1, l = 20, admix = 0.5, m1 = 0.05, m2 = 0.05,
          ne = 5000, p1 = 0.5, p2 = 0.5, nsims = "auto", tail_p = 0.2,
          keep_sims = FALSE)
  )
  hard <- suppressMessages(
    adsim(ngens = 33, k = 1.1, l = 5, admix = 0.5, m1 = 0, m2 = 0, ne = 200,
          p1 = 0.5, p2 = 0.9, nsims = "auto", tail_p = 0.2, keep_sims = FALSE)
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

test_that("resolving the tails is what decides nsims by default", {
  # Shape alone would settle these in a few hundred runs; supporting a test at
  # tail_p needs far more, and that is the binding requirement.
  easy <- suppressMessages(
    adsim(ngens = 20, k = 1, l = 20, admix = 0.5, m1 = 0.05, m2 = 0.05,
          ne = 5000, p1 = 0.5, p2 = 0.5, nsims = "auto", max_nsims = 2e5,
          keep_sims = FALSE)
  )
  shape_only <- suppressMessages(
    adsim(ngens = 20, k = 1, l = 20, admix = 0.5, m1 = 0.05, m2 = 0.05,
          ne = 5000, p1 = 0.5, p2 = 0.5, nsims = "auto", tail_p = 0.2,
          keep_sims = FALSE)
  )
  expect_gt(easy$conditions$nsims, 20 * shape_only$conditions$nsims)
  expect_equal(easy$conditions$nsims, tail_depth(0.5, 0.001, 0.2)$required,
               tolerance = 0.2)
  expect_true(easy$summary$tails_resolved)
})

test_that("tail depth follows 1 / tail_p when no atom bounds it", {
  spread <- seq(0.2, 0.8, length.out = 1000)   # no runs at 0 or 1
  shallow <- tail_depth(spread, tail_p = 0.01, tail_rse = 0.2)
  deep <- tail_depth(spread, tail_p = 0.001, tail_rse = 0.2)
  expect_equal(shallow$target, 0.01)
  expect_equal(deep$target, 0.001)
  # cost grows as 1 / tail_p
  expect_equal(deep$required / shallow$required, 10, tolerance = 0.02)
  expect_equal(deep$required, ceiling(0.999 / (0.001 * 0.04)))
})

test_that("an atom at a boundary caps the depth demanded", {
  # Half the runs lost and half fixed: no upper-tail p-value below about a
  # half exists, so asking for 0.001 cannot be met and must not be chased.
  split <- c(rep(0, 500), rep(1, 500))
  capped <- tail_depth(split, tail_p = 0.001, tail_rse = 0.2)
  expect_equal(capped$target, 0.5)
  expect_lt(capped$required, 100)

  # an atom on one side only still leaves the other tail to resolve
  one_side <- c(rep(0, 500), seq(0.1, 0.9, length.out = 500))
  single <- tail_depth(one_side, tail_p = 0.001, tail_rse = 0.2)
  expect_equal(single$target, 0.001)
})

test_that("auto samples deep enough to support a tail test", {
  res <- suppressMessages(
    adsim(ngens = 30, k = 1, l = 30, admix = 0.5, m1 = 0, m2 = 0, ne = 300,
          p1 = 0.5, p2 = 0.5, nsims = "auto", fixed_p0 = 0.5, tail_p = 0.01,
          max_nsims = 2e5, keep_sims = FALSE)
  )
  expect_true(res$summary$tails_resolved)
  expect_true(res$summary$converged)
  expect_equal(res$summary$tail_target, 0.01)
  expect_gte(res$conditions$nsims, ceiling(0.99 / (0.01 * 0.04)))
  expect_lte(res$summary$p_resolution, 0.01)
})

test_that("an absorbed scenario is not sampled past what it can resolve", {
  shallow <- suppressMessages(
    adsim(ngens = 1500, k = 1, l = 1500, admix = 0.5, m1 = 0, m2 = 0, ne = 50,
          p1 = 0.5, p2 = 0.5, nsims = "auto", precision = 0.02,
          tail_p = 0.001, fixed_p0 = 0.5, max_nsims = 2e5, keep_sims = FALSE)
  )
  # both boundaries carry atoms, so the target rises far above tail_p and the
  # run is not driven to the 250000 that 0.001 would otherwise demand
  expect_gt(shallow$summary$tail_target, 0.1)
  expect_lt(shallow$conditions$nsims, 50000)
  expect_true(shallow$summary$converged)
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

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

test_that("depth is set by where the observation falls, not by the scenario", {
  # The cost of resolving a tail probability depends on that probability alone,
  # so two very different nulls need the same depth for an observation at the
  # same percentile, and the same null needs very different depths for
  # observations at different percentiles.
  broad <- adsim_simulator(ngens = 33, k = 1.1, l = 5, admix = 0.5, m1 = 0,
                           m2 = 0, ne = 200, p1 = 0.5, p2 = 0.9,
                           nsims = 20000)$pF
  tight <- adsim_simulator(ngens = 20, k = 1, l = 20, admix = 0.5, m1 = 0.05,
                           m2 = 0.05, ne = 5000, p1 = 0.5, p2 = 0.5,
                           nsims = 20000)$pF

  at_decile <- function(x) {
    tail_depth(x, observed = stats::quantile(x, 0.9, names = FALSE))$required
  }
  expect_equal(at_decile(broad), at_decile(tight), tolerance = 0.1)
  expect_lt(at_decile(broad), 400)

  deep <- tail_depth(broad, observed = stats::quantile(broad, 0.999,
                                                       names = FALSE))
  expect_gt(deep$required, 50 * at_decile(broad))
})

test_that("no observation means no tail requirement", {
  spread <- seq(0.2, 0.8, length.out = 1000)
  none <- tail_depth(spread, observed = NULL)
  expect_equal(none$required, 0L)
  expect_true(is.na(none$target))
})

test_that("tail_p caps the cost of an extreme observation", {
  spread <- seq(0.2, 0.8, length.out = 1000)
  far <- tail_depth(spread, observed = 1, tail_p = 0.001, tail_rse = 0.2)
  expect_equal(far$target, 0.001)
  expect_equal(far$required, ceiling(0.999 / (0.001 * 0.04)))
  # a looser cap costs proportionally less
  loose <- tail_depth(spread, observed = 1, tail_p = 0.01, tail_rse = 0.2)
  expect_equal(loose$required, ceiling(0.99 / (0.01 * 0.04)))
})

test_that("an observation on an atom needs almost no depth", {
  # Half the runs fixed: an observation of 1 has a tail of about a half, which
  # is already the smallest p-value available in that direction.
  split <- c(rep(0, 500), rep(1, 500))
  on_atom <- tail_depth(split, observed = 1, tail_p = 0.001)
  expect_gt(on_atom$target, 0.4)
  expect_lt(on_atom$required, 100)
})

test_that("auto samples deep enough for the observation it is given", {
  res <- suppressMessages(
    adsim(ngens = 30, k = 1, l = 30, admix = 0.5, m1 = 0, m2 = 0, ne = 300,
          p1 = 0.5, p2 = 0.5, nsims = "auto", fixed_p0 = 0.5,
          observed = 0.56, tail_p = 0.001, max_nsims = 2e5, keep_sims = FALSE)
  )
  expect_true(res$summary$tails_resolved)
  expect_true(res$summary$converged)
  expect_false(is.na(res$summary$tail_target))
  # enough runs to pin the tail at that observation to about a fifth
  expect_gte(res$conditions$nsims,
             ceiling((1 - res$summary$tail_target) /
                       (res$summary$tail_target * 0.04)) * 0.9)
})

test_that("without an observation the shape criterion decides", {
  easy <- suppressMessages(
    adsim(ngens = 20, k = 1, l = 20, admix = 0.5, m1 = 0.05, m2 = 0.05,
          ne = 5000, p1 = 0.5, p2 = 0.5, nsims = "auto", keep_sims = FALSE)
  )
  hard <- suppressMessages(
    adsim(ngens = 33, k = 1.1, l = 5, admix = 0.5, m1 = 0, m2 = 0, ne = 200,
          p1 = 0.5, p2 = 0.9, nsims = "auto", keep_sims = FALSE)
  )
  expect_lt(easy$conditions$nsims, hard$conditions$nsims)
  expect_true(is.na(easy$summary$tail_target))
  expect_true(easy$summary$tails_resolved)
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

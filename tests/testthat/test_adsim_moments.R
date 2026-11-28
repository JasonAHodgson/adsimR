test_that("moments reduce to the classical drift decay", {
  # No migration and constant census size: var = p0 (1 - p0) [1 - (1 - 1/2N)^t]
  ne <- 250
  ngens <- 60
  got <- adsim_moments(ngens = ngens, k = 1, l = ngens, admix = 0.5,
                       m1 = 0, m2 = 0, ne = ne, p1 = 0.4, p2 = 0.4,
                       fixed_p0 = 0.4)
  expect_equal(got$mean, 0.4)
  expect_equal(got$var, 0.4 * 0.6 * (1 - (1 - (1 / (2 * ne)))^ngens))
})

test_that("moments match the simulator", {
  scenarios <- list(
    list(ngens = 50, k = 1, l = 50, admix = 0.5, m1 = 0, m2 = 0,
         ne = 200, p1 = 0.3, p2 = 0.9),
    list(ngens = 33, k = 1.1, l = 5, admix = 0.5, m1 = 0.01, m2 = 0.01,
         ne = 200, p1 = 0.5, p2 = 0.9),
    list(ngens = 40, k = 1, l = 40, admix = 0.37, m1 = 0, m2 = 0,
         ne = 101, p1 = 0.1, p2 = 0.95)
  )
  for (scenario in scenarios) {
    sims <- do.call(adsim_simulator, c(scenario, list(nsims = 40000)))
    expected <- do.call(adsim_moments, scenario)
    expect_equal(mean(sims$pF), expected$mean, tolerance = 0.02)
    expect_equal(stats::sd(sims$pF), expected$sd, tolerance = 0.05)
  }
})

test_that("census summaries describe a growing population", {
  got <- adsim_moments(ngens = 33, k = 1.1, l = 5, admix = 0.5,
                       m1 = 0, m2 = 0, ne = 200, p1 = 0.5, p2 = 0.9)
  # Growth begins after generation l, so the final size is the largest and the
  # harmonic mean lies well below it: drift is dominated by the early, small
  # generations.
  expect_equal(got$n_final, got$n_max)
  expect_gt(got$n_final, 200)
  expect_lt(got$n_harmonic, got$n_final)
  expect_gt(got$n_harmonic, 200)

  # A finite ceiling caps the trajectory.
  capped <- adsim_moments(ngens = 33, k = 1.1, l = 5, admix = 0.5,
                          m1 = 0, m2 = 0, ne = 200, p1 = 0.5, p2 = 0.9,
                          max_n = 500)
  expect_equal(capped$n_max, 500)
})

test_that("migration pulls the expectation toward the source", {
  # Allele absent from the admixed population, fixed in the second source.
  got <- adsim_moments(ngens = 20, k = 1, l = 20, admix = 1, m1 = 0, m2 = 0.1,
                       ne = 200, p1 = 0, p2 = 1)
  target <- 0
  for (i in seq_len(20)) target <- (0.9 * target) + 0.1
  expect_equal(got$mean, target)
  # The admixture expectation alone would give 0 here.
  expect_gt(got$mean, 0.8)
})

test_that("invalid conditions are rejected", {
  expect_error(adsim_moments(ngens = 10, k = 1, l = 20, admix = 0.5, m1 = 0,
                             m2 = 0, ne = 100, p1 = 0.5, p2 = 0.5))
  expect_error(adsim_moments(ngens = 10, k = 1, l = 5, admix = 0.5, m1 = 0.6,
                             m2 = 0.6, ne = 100, p1 = 0.5, p2 = 0.5))
})

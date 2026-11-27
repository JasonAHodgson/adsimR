test_that("drift variance matches the theoretical decay", {
  # Starting from a fixed p0 with no migration, the variance of allele
  # frequency after t generations is p0 (1 - p0) [1 - (1 - 1 / 2N)^t]. This
  # exercises the whole drift engine, including the removal of absorbed runs.
  for (params in list(c(ne = 100, ngens = 50), c(ne = 100, ngens = 400),
                      c(ne = 1000, ngens = 200))) {
    ne <- params[["ne"]]
    ngens <- params[["ngens"]]
    sims <- adsim_simulator(ngens = ngens, k = 1, l = ngens, admix = 0.5,
                            m1 = 0, m2 = 0, ne = ne, p1 = 0.5, p2 = 0.5,
                            nsims = 10000, fixed_p0 = 0.5)
    theory <- 0.25 * (1 - (1 - (1 / (2 * ne)))^ngens)
    expect_equal(stats::var(sims$pF), theory, tolerance = 0.06)
  }
})

test_that("effective population size is not capped by default", {
  # max_n once defaulted to 5000, which silently bounded every larger run.
  sims <- adsim_simulator(ngens = 200, k = 1, l = 200, admix = 0.5,
                          m1 = 0, m2 = 0, ne = 1e6, p1 = 0.5, p2 = 0.5,
                          nsims = 2000, fixed_p0 = 0.5)
  expect_equal(nrow(sims), 2000)
  expect_true(all(sims$pF >= 0 & sims$pF <= 1))
  # Drift is slow in a large population, so frequencies stay near p0.
  expect_equal(stats::var(sims$pF), 0.25 * (1 - (1 - (1 / 2e6))^200),
               tolerance = 0.1)

  # A finite ceiling must still be respected, and cannot sit below ne.
  expect_error(
    adsim_simulator(ngens = 10, k = 2, l = 1, admix = 0.5, m1 = 0, m2 = 0,
                    ne = 1000, p1 = 0.5, p2 = 0.5, nsims = 10, max_n = 500)
  )
})

test_that("absorbing states are preserved", {
  # With no migration an allele that is absent cannot appear, and one that is
  # fixed cannot be lost; dropping such runs from the sampler must not change
  # the values they report.
  absent <- adsim_simulator(ngens = 50, k = 1, l = 50, admix = 0.5,
                            m1 = 0, m2 = 0, ne = 100, p1 = 0, p2 = 0,
                            nsims = 200)
  expect_true(all(absent$pF == 0))

  fixed <- adsim_simulator(ngens = 50, k = 1, l = 50, admix = 0.5,
                           m1 = 0, m2 = 0, ne = 100, p1 = 1, p2 = 1,
                           nsims = 200)
  expect_true(all(fixed$pF == 1))

  # Migration from a polymorphic source makes neither state absorbing, so runs
  # must stay in the sampler and can move away from 0 and 1.
  rescued <- adsim_simulator(ngens = 50, k = 1, l = 50, admix = 0.5,
                             m1 = 0, m2 = 0.1, ne = 100, p1 = 0, p2 = 0.5,
                             nsims = 200, fixed_p0 = 0)
  expect_true(mean(rescued$pF) > 0)
})

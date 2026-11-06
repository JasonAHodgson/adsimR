test_that("p_value computes accurately", {
  test <- adsim(
    ngens = 10, k = 2, l = 5, admix = c(0.2, 0.5),
    m1 = 0.01, m2 = 0.01, ne = 100, p1 = 0.5,
    p2 = 0.9, nsims = 10
  )
  obs_freq <- 0.9
  p_out <- p_values(test, obs_freq)

  expect_equal(colnames(p_out), c("p1", "p2", "admix", "exp_freq", "z", "p"))
  # check that there is one row per unique combination of simulated parameters
  expect_equal(nrow(p_out), nrow(sim_overview(test)))
  # TODO this is not robust, needs comparison to another computation
  expect_equal(
    p_out$exp_freq,
    (p_out$admix * p_out$p1) + ((1 - p_out$admix) * p_out$p2)
  )
})

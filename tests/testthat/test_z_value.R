test_that("z_value computes accurately", {
  test <- adsim(
    ngens = 10, k = 2, l = 5, admix = c(0.2, 0.5),
    m1 = 0.01, m2 = 0.01, ne = 100, p1 = 0.5,
    p2 = 0.9, nsims = 10
  )
  obs_freq <- 0.9
  z_output <- z_value(test, obs_freq)

  expect_equal(colnames(z_output), c("p1", "p2", "admix", "exp_freq", "z"))
  # check that there is one row per unique combination of simulated parameters
  expect_equal(nrow(z_output), nrow(sim_overview(test)))
  # TODO this is not robust, needs comparison to another computation
  expect_equal(
    z_output$exp_freq,
    (z_output$admix * z_output$p1) + ((1 - z_output$admix) * z_output$p2)
  )
})

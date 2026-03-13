n_reps <- 10

test_that("test output dimensions", {
  admix_props <- c(0.2, 0.5)
  test <- adsim(ngens = 10, k = 2, l = 5, admix = admix_props,
                m1 = 0.01, m2 = 0.01, ne = 100, p1 = 0.5,
                p2 = 0.9, nsims = n_reps)
  expect_equal(nrow(test), as.integer(length(admix_props)*n_reps))
})



test_that("adsim output is compatible with adsim_simulator output", {
t_1 <- adsim_simulator(ngens = 10, k = 2, l = 5, admix = 0.2,
               m1 = 0.01, m2 = 0.01, ne = 100, p1 = 0.5,
               p2 = 0.9, nsims = n_reps)
t_2 <- adsim_simulator(ngens = 10, k = 2, l = 5, admix = 0.4,
                         m1 = 0.01, m2 = 0.01, ne = 100, p1 = 0.5,
                         p2 = 0.9, nsims = n_reps)
t_3 <- adsim_simulator(ngens = 10, k = 2, l = 5, admix = 0.6,
                         m1 = 0.01, m2 = 0.01, ne = 100, p1 = 0.5,
                         p2 = 0.9, nsims = n_reps)

adsim_simulator_output <- bind_rows(t_1, t_2, t_3)

adsim_output <- adsim(ngens = 10, k = 2, l = 5, admix = c(0.2, 0.4, 0.6),
             m1 = 0.01, m2 = 0.01, ne = 100, p1 = 0.5,
             p2 = 0.9, nsims = n_reps)
expect_equal(adsim_output$admix, adsim_simulator_output$admix)
})

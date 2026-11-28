n_reps <- 10

test_that("test output dimensions", {
  admix_props <- c(0.2, 0.5)
  test <- adsim(ngens = 10, k = 2, l = 5, admix = admix_props,
                m1 = 0.01, m2 = 0.01, ne = 100, p1 = 0.5,
                p2 = 0.9, nsims = n_reps)
  expect_equal(nrow(test$sims), as.integer(length(admix_props)*n_reps))
  expect_equal(nrow(test$conditions), length(admix_props))
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

adsim_output <- as.data.frame(adsim(ngens = 10, k = 2, l = 5,
             admix = c(0.2, 0.4, 0.6),
             m1 = 0.01, m2 = 0.01, ne = 100, p1 = 0.5,
             p2 = 0.9, nsims = n_reps))
expect_equal(adsim_output$admix, adsim_simulator_output$admix)
})


test_that("migration moves the allele frequency toward the source", {
  # The allele is absent from the admixed population (admix = 1 with p1 = 0)
  # and fixed in the second parent population, so any copy observed at the end
  # must have arrived by migration.
  no_mig <- as.data.frame(adsim(ngens = 20, k = 1.05, l = 5, admix = 1,
                  m1 = 0, m2 = 0, ne = 200, p1 = 0, p2 = 1, nsims = 200))
  expect_equal(mean(no_mig$pF), 0)

  with_mig <- as.data.frame(adsim(ngens = 20, k = 1.05, l = 5, admix = 1,
                    m1 = 0, m2 = 0.1, ne = 200, p1 = 0, p2 = 1, nsims = 200))

  # Binomial sampling is unbiased and the migration step is linear in p, so the
  # mean final frequency tracks the deterministic continent-island recursion,
  # applied once per simulated generation.
  expected <- 0
  for (i in seq_len(20)) expected <- (0.9 * expected) + (0.1 * 1)
  expect_equal(mean(with_mig$pF), expected, tolerance = 0.03)
})

test_that("migration rates outside the unit interval are rejected", {
  # adsim skips combinations it cannot simulate and only fails when none
  # remain, so a single impossible combination is an error.
  expect_error(
    adsim(ngens = 5, k = 1, l = 1, admix = 0.5, m1 = 0.6, m2 = 0.6,
          ne = 50, p1 = 0.5, p2 = 0.5, nsims = 1),
    "No valid combinations"
  )
  expect_error(
    adsim(ngens = 5, k = 1, l = 1, admix = 0.5, m1 = 1.5, m2 = 0,
          ne = 50, p1 = 0.5, p2 = 0.5, nsims = 1),
    "No valid combinations"
  )

  # When only some combinations are impossible the rest still run.
  expect_warning(
    kept <- suppressMessages(as.data.frame(
      adsim(ngens = 5, k = 1, l = 1, admix = 0.5, m1 = c(0.1, 0.6),
            m2 = 0.6, ne = 50, p1 = 0.5, p2 = 0.5, nsims = 10)
    )),
    "combination\\(s\\) skipped"
  )
  expect_equal(unique(kept$m1), 0.1)
})

# TODO make parameter checks have expected error and warning messages
test_that("adsim warnings",{

  expect_error(adsim(ngens = 10, k = 2, l = 10, admix = admix_list,
                     m1 = 0.01, m2 = 0.01, ne = 100, p1 = 0.5,
                     p2 = 0.9, nsims = 10),

  )

})

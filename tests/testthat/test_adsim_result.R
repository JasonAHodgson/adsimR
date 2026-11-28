make_result <- function(...) {
  suppressMessages(
    adsim(ngens = 25, k = 1.1, l = 5, admix = c(0.3, 0.5), m1 = 0, m2 = 0.01,
          ne = 150, p1 = 0.4, p2 = 0.9, nsims = 500, ...)
  )
}

test_that("adsim returns a three part result object", {
  res <- make_result()
  expect_s3_class(res, "adsim_result")
  expect_named(res, c("conditions", "sims", "summary"))

  expect_equal(nrow(res$conditions), 2)
  expect_equal(nrow(res$summary), 2)
  expect_equal(nrow(res$sims), 1000)
  expect_named(res$sims, c("scenario", "sim", "p0", "pF"))

  # conditions are recorded once per scenario, not once per run
  expect_equal(res$conditions$admix, c(0.3, 0.5))
  expect_true(all(c("n_final", "n_max", "n_harmonic") %in%
                    names(res$conditions)))
  expect_false(is.null(attr(res, "created")))
  expect_false(is.null(attr(res, "call")))
})

test_that("summary statistics describe the stored runs", {
  res <- make_result()
  for (i in seq_len(nrow(res$summary))) {
    runs <- res$sims[res$sims$scenario == i, ]
    row <- res$summary[i, ]
    expect_equal(row$mean_pF, mean(runs$pF))
    expect_equal(row$var_pF, stats::var(runs$pF))
    expect_equal(row$sd_pF, stats::sd(runs$pF))
    expect_equal(row$prop_fixed, mean(runs$pF >= 1))
    expect_equal(row$prop_lost, mean(runs$pF <= 0))
    expect_equal(row$median, stats::median(runs$pF))
    expect_equal(row$se_mean_pF, stats::sd(runs$pF) / sqrt(nrow(runs)))
  }
  # the three absorption proportions partition the runs
  expect_equal(res$summary$prop_fixed + res$summary$prop_lost +
                 res$summary$prop_segregating, rep(1, nrow(res$summary)))
})

test_that("fst equals the loss of heterozygosity under drift from a fixed p0", {
  # The identity fst == 1 - het_end / het_start needs the mean to be unchanged
  # (no migration) and the starting frequency to have no variance across runs.
  # A sampled founding generation gives var(p0) > 0, so fixed_p0 is required
  # for the identity to be exact. This is why the object reports the two
  # heterozygosities raw rather than as a proportion retained.
  drift <- suppressMessages(
    adsim(ngens = 100, k = 1, l = 100, admix = 0.5, m1 = 0, m2 = 0,
          ne = 50, p1 = 0.5, p2 = 0.5, nsims = 2000, fixed_p0 = 0.5,
          keep_sims = FALSE)
  )
  expect_equal(drift$summary$sd_p0, 0)
  # var_pF is the sample variance (n - 1) while het_end is a plain mean over
  # runs (n), so the two sides differ by the Bessel correction, of relative
  # size 1 / nsims.
  expect_equal(drift$summary$fst,
               1 - (drift$summary$het_end / drift$summary$het_start),
               tolerance = 2e-3)
  expect_gt(drift$summary$fst, 0)
  expect_lt(drift$summary$fst, 1)

  # Run far past 4N and differentiation is essentially complete.
  done <- suppressMessages(
    adsim(ngens = 1500, k = 1, l = 1500, admix = 0.5, m1 = 0, m2 = 0,
          ne = 50, p1 = 0.5, p2 = 0.5, nsims = 2000, fixed_p0 = 0.5,
          keep_sims = FALSE)
  )
  expect_lt(done$summary$prop_segregating, 0.005)
  expect_equal(done$summary$fst, 1, tolerance = 0.02)
})

test_that("the analytic cross-check agrees with the simulation", {
  res <- make_result()
  expect_equal(res$summary$sd_ratio, rep(1, nrow(res$summary)),
               tolerance = 0.15)
  expect_equal(res$summary$mean_expected,
               sapply(seq_len(nrow(res$conditions)), function(i) {
                 cond <- res$conditions[i, ]
                 adsim_moments(cond$ngens, cond$k, cond$l, cond$admix, cond$m1,
                               cond$m2, cond$ne, cond$p1, cond$p2,
                               max_n = cond$max_n)$mean
               }))
})

test_that("keep_sims = FALSE drops the runs but keeps the summary", {
  res <- make_result(keep_sims = FALSE)
  expect_null(res$sims)
  expect_equal(nrow(res$summary), 2)
  expect_equal(res$summary$nsims, c(500, 500))
  expect_error(as.data.frame(res), "keep_sims = TRUE")
  expect_lt(utils::object.size(res), utils::object.size(make_result()))
})

test_that("as.data.frame rebuilds the flat one row per run format", {
  res <- make_result()
  flat <- as.data.frame(res)
  expect_equal(nrow(flat), 1000)
  expect_true(all(c("scenario", "sim", "p0", "pF", "ngens", "k", "l", "admix",
                    "m1", "m2", "ne", "p1", "p2") %in% names(flat)))
  # every run carries the conditions of its own scenario
  first <- flat[flat$scenario == 1, ]
  expect_equal(unique(first$admix), 0.3)
  expect_equal(sort(flat$pF), sort(res$sims$pF))
})

test_that("downstream functions accept the result object", {
  res <- make_result()
  from_object <- sim_overview(res)
  from_frame <- sim_overview(as.data.frame(res))
  expect_equal(from_object, from_frame)
  expect_equal(nrow(from_object), 2)

  # z_value and p_values reach sim_overview, so they accept it too
  expect_equal(nrow(p_values(res, 0.8)), 2)
  expect_equal(nrow(z_value(res, 0.8)), 2)
})

test_that("print and summary methods behave", {
  res <- make_result()
  expect_identical(summary(res), res$summary)
  expect_output(print(res), "<adsim_result>")
  expect_output(print(res), "2 scenario")
})

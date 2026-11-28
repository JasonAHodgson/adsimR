drift_sims <- function(nsims = 4000, ...) {
  suppressMessages(
    adsim(ngens = 30, k = 1, l = 30, admix = 0.5, m1 = 0, m2 = 0, ne = 200,
          p1 = 0.5, p2 = 0.5, nsims = nsims, fixed_p0 = 0.5, ...)
  )
}

test_that("the p-value counts runs beyond the observation", {
  sims <- drift_sims()
  got <- adsim_pvalue(sims, observed = 0.6, alternative = "greater")
  runs <- sims$sims$pF
  expect_equal(got$n_beyond, sum(runs >= 0.6))
  expect_equal(got$p_value, (sum(runs >= 0.6) + 1) / (length(runs) + 1))
  expect_equal(got$nsims, length(runs))
})

test_that("ties count as extreme", {
  # A run exactly at the observed value is as extreme as the observation, and
  # dropping such runs would be anti-conservative where the atoms at 0 and 1
  # sit.
  sims <- drift_sims()
  sims$sims$pF <- c(rep(0.5, 100), rep(0.9, 100))
  sims$sims <- sims$sims[seq_len(200), ]
  sims$sims$scenario <- 1
  sims$conditions <- sims$conditions[1, ]
  got <- adsim_pvalue(sims, observed = 0.5, alternative = "greater")
  expect_equal(got$n_beyond, 200)
})

test_that("a p-value of zero is never reported", {
  sims <- drift_sims()
  got <- adsim_pvalue(sims, observed = 1, alternative = "greater")
  expect_gt(got$p_value, 0)
  expect_equal(got$p_value, 1 / (got$nsims + 1))
  expect_equal(got$p_resolution, 1 / (got$nsims + 1))
  expect_false(got$resolved)
})

test_that("two sided is the default and doubles the smaller tail", {
  sims <- drift_sims()
  both <- adsim_pvalue(sims, observed = 0.62)
  expect_equal(both$alternative, "two.sided")

  greater <- adsim_pvalue(sims, observed = 0.62, alternative = "greater")
  less <- adsim_pvalue(sims, observed = 0.62, alternative = "less")
  expect_equal(both$p_value, min(1, 2 * min(greater$p_value, less$p_value)))
  expect_lte(both$p_value, 1)
})

test_that("direction behaves as it should", {
  sims <- drift_sims()
  high <- adsim_pvalue(sims, observed = 0.75, alternative = "greater")
  low <- adsim_pvalue(sims, observed = 0.75, alternative = "less")
  # a high observation is unusual from above and unremarkable from below
  expect_lt(high$p_value, low$p_value)
  expect_gt(low$p_value, 0.5)
})

test_that("the confidence interval brackets the p-value", {
  sims <- drift_sims()
  got <- adsim_pvalue(sims, observed = 0.58)
  expect_lte(got$p_lower, got$p_value)
  expect_gte(got$p_upper, got$p_value)
  expect_gte(got$p_lower, 0)
  expect_lte(got$p_upper, 1)
})

test_that("one row is returned per scenario", {
  sims <- suppressMessages(
    adsim(ngens = 30, k = 1, l = 30, admix = c(0.3, 0.7), m1 = 0, m2 = 0,
          ne = 200, p1 = 0.2, p2 = 0.9, nsims = 2000)
  )
  got <- adsim_pvalue(sims, observed = 0.8)
  expect_equal(nrow(got), 2)
  expect_equal(got$scenario, c(1, 2))
  # the same observation is judged against each scenario's own null
  expect_false(isTRUE(all.equal(got$p_value[1], got$p_value[2])))
})

test_that("the test is valid when observations come from the null", {
  # Drawing observations from the null itself, a valid test rejects at no more
  # than the nominal rate.
  sims <- drift_sims(nsims = 20000)
  null <- sims$sims$pF
  draws <- sample(null, 2000)
  p <- vapply(draws, function(o) {
    upper <- sum(null >= o)
    lower <- sum(null <= o)
    min(1, 2 * ((min(upper, lower) + 1) / (length(null) + 1)))
  }, numeric(1))
  expect_lte(mean(p <= 0.05), 0.06)
  expect_lte(mean(p <= 0.01), 0.02)
})

test_that("inputs are checked", {
  sims <- drift_sims(nsims = 500)
  expect_error(adsim_pvalue(sims$summary, observed = 0.5), "adsim_result")
  expect_error(adsim_pvalue(sims, observed = 1.5))
  expect_error(adsim_pvalue(sims, observed = NA_real_))
  lean <- drift_sims(nsims = 500, keep_sims = FALSE)
  expect_error(adsim_pvalue(lean, observed = 0.5), "keep_sims = TRUE")
})

test_that("sim_overview summarises correctly",{

  admix_list <- c(0.2, 0.5)
  gens_list <- c(10, 20)
  test <- adsim(ngens = gens_list, k = 2, l = 5, admix = admix_list,
  m1 = 0.01, m2 = 0.01, ne = 100, p1 = 0.5,
  p2 = 0.9, nsims = 10)

  overview_obj <- sim_overview(test)

  # check the number of rows are equal to the number of unique
  # parameter combinations (length of admix_list * length of gens_list)
  expect_equal(nrow(overview_obj), length(admix_list) * length(gens_list))

  # check that the columns are as expected
  expected_cols <- c("ngens", "k", "l", "admix", "m1", "m2", "ne", "p1", "p2",
                     "mean_pF", "sd_pF", "Q1_pF", "Q2_pF",
                     "Q3_pF", "Q4_pF", "Q5_pF", "Q6_pF",
                     "param_group")
  expect_equal(colnames(overview_obj), expected_cols)

  # check that subsets are correctly calculated
  for (i in 1:nrow(overview_obj)) {
    subset_df <- test[test$ngens == overview_obj$ngens[i] &
                       test$admix == overview_obj$admix[i], ]
    expected_mean <- mean(subset_df$pF)
    expect_equal(overview_obj$mean_pF[i], expected_mean)
    expected_sd <- stats::sd(subset_df$pF)
    expect_equal(overview_obj$sd_pF[i], expected_sd)
  }
})

test_that("sim_overview when every parameter is > length 1",{
  admix_list <- c(0.2, 0.5, 0.7)
  gens_list <- c(10, 20)
  k_list <- c(1, 2)
  l_list <- c(5, 7)
  m1_list <- c(0.01, 0.02)
  m2_list <- c(0.01, 0.02)
  ne_list <- c(100, 200, 300)
  p1_list <- c(0.5, 0.6)
  p2_list <- c(0.9, 0.8)
  n_sims <- 4

  test <- adsim(ngens = gens_list, k = k_list, l = l_list, admix = admix_list,
                 m1 = m1_list, m2 = m2_list, ne = ne_list, p1 = p1_list,
                 p2 = p2_list, nsims = n_sims)

  # check the number of rows are equal to the number of unique parameter combinations
  expected_rows <- length(admix_list) * length(gens_list) * length(k_list) *
    length(l_list) * length(m1_list) * length(m2_list) * length(ne_list) *
    length(p1_list) * length(p2_list)
  overview_obj <- sim_overview(test)
  expect_equal(nrow(overview_obj), expected_rows)

})

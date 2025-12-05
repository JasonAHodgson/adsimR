#' sim_overview
#'
#' This function creates an overview summary of adsim output
#' @param adsim_df  A data frame containing three columns:
#' - Sim: a numeric index containing the simulation run
#' - p0: a numeric containing the starting allele following the first
#'   generation of admixture.
#' - pF: a numeric containing the final allele frequency of each simulation.
#' @return a data frame summarizing the mean, standard deviation, and percentiles
#' @export
#' @examples
#' test <- adsim(ngens = 10, k = 2, l = 5, admix = c(0.2, 0.5),
#' m1 = 0.01, m2 = 0.01, ne = 100, p1 = 0.5,
#' p2 = 0.9, nsims = 10)
#'
#' sim_overview(test)

sim_overview <- function(adsim_df){
  sum_adsim_df <- adsim_df %>%
    dplyr::group_by(ngens, k, l, admix, m1, m2, ne, p1, p2) %>%
    dplyr::summarise(
      mean_pF = mean(pF),
      sd_pF   = stats::sd(pF),
      # calculate percentiles
      Q1_pF = stats::quantile(pF, 0.05),
      Q2_pF = stats::quantile(pF, 0.10),
      Q3_pF = stats::quantile(pF, 0.25),
      Q4_pF = stats::quantile(pF, 0.75),
      Q5_pF = stats::quantile(pF, 0.90),
      Q6_pF = stats::quantile(pF, 0.95),
      .groups = "drop"
    )
  # Create vector P1..PN where N = number of rows
    sum_adsim_df$param_group <- paste0("P", seq_len(nrow(sum_adsim_df)))
    return(sum_adsim_df)
}

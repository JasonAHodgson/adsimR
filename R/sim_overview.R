#' sim_overview
#'
#' This function creates an overview summary of `adsim` output.
#'
#' @param adsim_df  A data frame containing the output of an `adsim` simulation,
#'   with the following columns:
#' - Sim: a numeric index containing the simulation run
#' - p0: a numeric containing the starting allele frequency following the first
#'   generation of admixture.
#' - pF: a numeric containing the final allele frequency of each simulation.
#' - ngens: Number of generations simualted.
#' - k: population growth parameter
#' - l: number of generations prior to population growth.
#' - admix: Proportion of ancestry contributed by the first parent
#' population.
#' - ne: effective population size of generation 0.
#' - p1: Initial allele frequency in the first parent population.
#' - p2: Initial allele frequency in the second parent population.
#' - m1: Migration rate from the first parent population.
#' - m2: Migration rate from the second parent population.
#' @return A data frame summarizing, for each unique set of parameters, the
#'   mean, standard deviation, and percentiles across the set of n simulations.
#'   Columns of the returned data frame are:
#'   - ngens: Number of generations simualted.
#'   - k: population growth parameter
#'   - l: number of generations prior to population growth.
#'   - admix: Proportion of ancestry contributed by the first parent
#'   population.
#'   - m1: Migration rate from the first parent population.
#'   - m2: Migration rate from the second parent population.
#'   - ne: effective population size of generation 0.
#'   - p1: Initial allele frequency in the first parent population.
#'   - p2: Initial allele frequency in the second parent population.
#'   - mean_pF: Mean final allele frequency across simulations.
#'   - sd_pF: Standard deviation of final allele frequency across simulations.
#'   - Q1_pF to Q6_pF: Percentiles of final allele frequency across simulations.
#'   - param_group: A unique identifier for each combination of parameters,
#'     i.e. each unique set of simulations.
#' @export
#' @examples
#' test <- adsim(ngens = 10, k = 2, l = 5, admix = c(0.2, 0.5),
#' m1 = 0.01, m2 = 0.01, ne = 100, p1 = 0.5,
#' p2 = 0.9, nsims = 10)
#'
#' sim_overview(test)

sim_overview <- function(adsim_df){
  sum_adsim_df <- adsim_df %>%
    dplyr::group_by(.data$ngens, .data$k, .data$l, .data$admix,
                    .data$m1, .data$m2, .data$ne, .data$p1,
                    .data$p2) %>%
    dplyr::summarise(
      mean_pF = mean(.data$pF),
      sd_pF   = stats::sd(.data$pF),
      # calculate percentiles
      Q1_pF = stats::quantile(.data$pF, 0.05),
      Q2_pF = stats::quantile(.data$pF, 0.10),
      Q3_pF = stats::quantile(.data$pF, 0.25),
      Q4_pF = stats::quantile(.data$pF, 0.75),
      Q5_pF = stats::quantile(.data$pF, 0.90),
      Q6_pF = stats::quantile(.data$pF, 0.95),
      .groups = "drop"
    )
  # Create vector P1..PN where N = number of rows
    sum_adsim_df$param_group <- paste0("P", seq_len(nrow(sum_adsim_df)))
    return(sum_adsim_df)
}

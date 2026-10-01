#' p_value
#'
#' This function calculates p-values for observed vs expected allele
#' frequencies.
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
#'   population.
#' - ne: effective population size of generation 0.
#' - p1: Initial allele frequency in the first parent population.
#' - p2: Initial allele frequency in the second parent population.
#' - m1: Migration rate from the first parent population.
#' - m2: Migration rate from the second parent population.
#' @param observed_freq Observed allele frequency
#' @return A data frame returning a p-value and a z-score for the observed
#'   allele frequency compared to the expected allele frequency derived from the
#'   simulated dataset. The output contains one row for each unique set of
#'   simulated parameters, with the following columns:
#'   - p1: Initial allele frequency in the first parent population.
#'   - p2: Initial allele frequency in the second parent population.
#'   - admix: Proportion of ancestry contributed by the first parent
#'   population.
#'   - exp_freq: Expected allele frequency based on the admixture proportion and
#'   the initial allele frequencies of the two parent populations.
#'   - z: z-score calculated as (observed_freq - exp_freq) / sd_pF,
#'   where exp_freq is based on the admixture proportion and the initial allele
#'   frequencies of the two populations, and sd_pF is the standard deviation of
#'   the final allele frequencies across all simulations.
#'   - p: p-value calculated as 2*(1 - pnorm(abs(z))), where z is the z-score
#'   described in the documentation for [z_value()].
#' @export
#' @examples
#' test <- adsim(ngens = 10, k = 2, l = 5, admix = c(0.2, 0.5),
#' m1 = 0.01, m2 = 0.01, ne = 100, p1 = 0.5,
#' p2 = 0.01, nsims = 10)
#'
#' obs_freq <- 0.9
#'
#' p_values(test, obs_freq)

p_values <- function(adsim_df, observed_freq){

  standard <- sim_overview(adsim_df) %>%
    dplyr::mutate(obs_freq = observed_freq)

  standard <- standard %>%
    dplyr::group_by(.data$p1, .data$p2, .data$admix) %>%
    dplyr::transmute(
      #obs_freq = .data$mean_pF,
      exp_freq = (.data$admix * .data$p1) + ((1-.data$admix)* .data$p2),
      z = (.data$obs_freq - .data$exp_freq)/.data$sd_pF,
      p = 2*(1-stats::pnorm(abs(.data$z)))
    )
  return(standard)
}

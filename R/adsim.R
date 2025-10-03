#' adsim
#'
#' This function simulates evolution of an allele through genetic drift in an
#' admixed population. The user must specify the starting allele frequency in
#' both parent populations, the number of generations to simulate, the starting
#' effective population size, the admixture proportion, the population growth
#' parameters, the migration rate from both parent populations, and the number
#' of simulations to perform. The function returns a data frame containing the
#' final allele frequency of each simulation.The function accepts lists as input
#' as well.
#'
#' @param adsim_simulator Function that performs a single simulation run.
#' @param ngens Number of generations to simulate.
#' @param k population growth parameter.
#' @param l number of generations prior to population growth.
#' @param admix Proportion of ancestry contributed by the first parent population.
#' @param ne effective population size of generation 0.
#' @param p1 Initial allele frequency in the first parent population.
#' @param p2 Initial allele frequency in the second parent population.
#' @param m1 Migration rate from the first parent population.
#' @param m2 Migration rate from the second parent population.
#' @param nsims Number of simulations to perform.
#' @return A data frame containing three columns:
#' - Sim: a numeric index containing the simulation run
#' - p0: a numeric containing the starting allele following the first generation of admixture.
#' - pF: a numeric containing the final allele frequency of each simulation.
#' @export

adsim <- function(ngens, k, l, admix,
                      m1, m2, ne, p1, p2, nsims)
  {
  # Input validation for vectors
  stopifnot(all(is.numeric(ngens), ngens > 0))
  stopifnot(all(is.numeric(k), k > 0))
  stopifnot(all(is.numeric(l), l >= 0, l < ngens))
  stopifnot(all(is.numeric(admix), admix >= 0, admix <= 1))
  stopifnot(all(is.numeric(m1), m1 >= 0))
  stopifnot(all(is.numeric(m2), m2 >= 0))
  stopifnot(all(is.numeric(ne), ne > 0))
  stopifnot(all(is.numeric(p1), p1 >= 0, p1 <= 1))
  stopifnot(all(is.numeric(p2), p2 >= 0, p2 <= 1))
  stopifnot(all(is.numeric(nsims), nsims > 0))

  # Create parameter combinations
  param_combos <- expand.grid(
    ngens = ngens, k = k, l = l, admix = admix,
    m1 = m1, m2 = m2, ne = ne, p1 = p1, p2 = p2, nsims = nsims
  )

  # Sanity check - compute the number of combinations of parameters
  n_combinations <- length(ngens) * length(k) * length(l) * length(admix) *
    length(m1) * length(m2) * length(ne) * length(p1) * length(p2) * length(nsims)
  total_sims <- n_combinations * nsims


  # Function to run mapply and combine results
  combine_sims <- function(ngens, k, l, admix, m1, m2, ne, p1, p2, nsims) {
    results_list <- mapply(adsim_simulator, ngens, k, l, admix, m1, m2, ne, p1, p2, nsims, SIMPLIFY = FALSE)
    results <- do.call(rbind, results_list)
    return(results)
  }

  # Run simulations and combine results
  return(combine_sims(
    param_combos$ngens, param_combos$k, param_combos$l, param_combos$admix,
    param_combos$m1, param_combos$m2, param_combos$ne, param_combos$p1,
    param_combos$p2, param_combos$nsims
  ))
}


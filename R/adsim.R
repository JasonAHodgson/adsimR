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
#' @param ngens Number of generations to simulate.
#' @param k population growth parameter.
#' @param l number of generations prior to population growth.
#' @param admix Proportion of ancestry contributed by the first parent
#'   population.
#' @param ne effective population size of generation 0.
#' @param p1 Initial allele frequency in the first parent population.
#' @param p2 Initial allele frequency in the second parent population.
#' @param m1 Per-generation migration rate from the first parent population,
#'   that is the proportion of gene copies in the admixed population replaced
#'   each generation by copies drawn from that population. Must lie in
#'   \[0, 1\], and `m1 + m2` must not exceed 1.
#' @param m2 Per-generation migration rate from the second parent population,
#'   defined as for `m1`.
#' @param nsims Number of simulations to perform.
#' @return A data frame containing three columns:
#' - Sim: a numeric index containing the simulation run
#' - p0: a numeric containing the starting allele following the first generation
#'   of admixture.
#' - pF: a numeric containing the final allele frequency of each simulation.
#' @details Generation 0 is formed by sampling `admix * ne` individuals from
#'   the first parent population and the remainder from the second. Each
#'   subsequent generation is produced in two steps. Migration acts first on
#'   the gene pool, following the continent-island model, so that the frequency
#'   entering the drift step is
#'   \deqn{p^{*} = (1 - m_1 - m_2)\,p + m_1 p_1 + m_2 p_2,}
#'   and genetic drift then follows, with \eqn{2N} gene copies drawn
#'   binomially at frequency \eqn{p^{*}}. Migration therefore changes allele
#'   frequency but not census size, which is governed by `k` and `l` alone and
#'   is capped at 5000 individuals.
#' @export

adsim <- function(ngens, k, l, admix,
                      m1, m2, ne, p1, p2, nsims)
  {
  # Input validation for vectors
  stopifnot(all(is.numeric(ngens), ngens > 0))
  stopifnot(all(is.numeric(k), k > 0))
  stopifnot(all(is.numeric(l), l >= 0, l < ngens))
  stopifnot(all(is.numeric(admix), admix >= 0, admix <= 1))
  stopifnot(all(is.numeric(m1), m1 >= 0, m1 <= 1))
  stopifnot(all(is.numeric(m2), m2 >= 0, m2 <= 1))
  stopifnot(all(is.numeric(ne), ne > 0))
  stopifnot(all(is.numeric(p1), p1 >= 0, p1 <= 1))
  stopifnot(all(is.numeric(p2), p2 >= 0, p2 <= 1))
  stopifnot(all(is.numeric(nsims), nsims > 0))

  # Create parameter combinations
  param_combos <- expand.grid(
    ngens = ngens, k = k, l = l, admix = admix,
    m1 = m1, m2 = m2, ne = ne, p1 = p1, p2 = p2, nsims = nsims
  )

  # m1 and m2 are crossed, so the combined rate has to be checked per
  # combination rather than on the supplied vectors.
  if (any(param_combos$m1 + param_combos$m2 > 1)) {
    stop("m1 + m2 must not exceed 1: migration rates are the proportions of ",
         "the gene pool replaced each generation.")
  }

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


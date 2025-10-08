#' adsim_simulator
#'
#' This function simulates evolution of an allele through genetic drift in an
#' admixed population. The user must specify the starting allele frequency in
#' both parent populations, the number of generations to simulate, the starting
#' effective population size, the admixture proportion, the population growth
#' parameters, the migration rate from both parent populations, and the number
#' of simulations to perform. The function returns a data frame containing the
#' final allele frequency of each simulation.
#'
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

adsim_simulator <- function(
    ngens,      # generations since admixture
    k,          # population growth factor
    l,          # generations before growth
    admix,      # admixture proportion
    m1,         # migration rate from population 1
    m2,         # migration rate from population 2
    ne,         # effective population size at generation 0
    p1,         # allele frequency in population 1
    p2,         # allele frequency in population 2
    nsims       # number of simulations
) {
  max_n <- 5000

  # Input validation
  stopifnot(is.numeric(ngens), ngens > 0)
  stopifnot(is.numeric(k), k > 0)
  stopifnot(is.numeric(l), l >= 0, l < ngens)
  stopifnot(is.numeric(admix), admix >= 0, admix <= 1)
  stopifnot(is.numeric(m1), m1 >= 0)
  stopifnot(is.numeric(m2), m2 >= 0)
  stopifnot(is.numeric(ne), ne > 0)
  stopifnot(is.numeric(p1), p1 >= 0, p1 <= 1)
  stopifnot(is.numeric(p2), p2 >= 0, p2 <= 1)
  stopifnot(is.numeric(nsims), nsims > 0)

  choose_allele <- function(n, p) {
    sum(runif(2 * n) <= p)
  }

  results <- data.frame(Sim = integer(), p0 = numeric(), pF = numeric(),
                        pF = numeric(), ngens = integer(),
                        k = numeric(), l = integer(),
                        admix = numeric(), m1 = numeric(),
                        m2 = numeric(), ne = integer(),
                        p1 = numeric(), p2 = numeric(),
                        stringsAsFactors = FALSE)

    for (sim in 1:nsims) {
      raw_p1 <- admix * ne
      n_p1 <- if (raw_p1 %% 1 != 0) {
        int <- floor(raw_p1)
        fp <- raw_p1 - int
        if (runif(1) > fp) int else int + 1
      } else raw_p1
      n_p2 <- ne - n_p1

      P1_count <- choose_allele(n_p1, p1)
      P2_count <- choose_allele(n_p2, p2)
      A_initial <- (P1_count + P2_count) / (2 * ne)

      g_count <- 1
      old_n <- ne
      old_A <- A_initial

      while (g_count < ngens) {
        new_n <- if (g_count > l) floor(old_n * k) else ne
        new_n <- min(new_n, max_n)

        A_count <- choose_allele(new_n, old_A)

        if (m1 > 0) {
          A_count <- A_count + choose_allele(m1, p1)
          new_n <- new_n + m1
        }

        if (m2 > 0) {
          A_count <- A_count + choose_allele(m2, p2)
          new_n <- new_n + m2
        }

        old_A <- A_count / (2 * new_n)
        old_n <- new_n
        g_count <- g_count + 1

        if (g_count == ngens) {
          results <- rbind(results, data.frame(Sim = sim, p0 = A_initial, pF = old_A, ngens = ngens, k = k,
                                               l = l, admix = admix, m1 = m1, m2 = m2, ne = ne, p1 = p1,
                                               p2 = p2, stringsAsFactors = FALSE) )
        }
      }

  }


  return(results)
}

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
#' @param ngens Number of generations to simulate after admixture.
#' @param k Population growth multiplier per generation.
#' @param l Number of generations prior to population growth.
#' @param admix Proportion of ancestry contributed by the first parent population.
#' @param m1 Migration proportion from the first parent population (0 to 1).
#' @param m2 Migration proportion from the second parent population (0 to 1).
#' @param ne Effective population size of generation 0.
#' @param p1 Initial allele frequency in the first parent population.
#' @param p2 Initial allele frequency in the second parent population.
#' @param nsims Number of simulations to perform.
#' @param fixed_p0 Optional fixed starting allele frequency for every run.
#' @param max_n Maximum simulated population size. Default 5000.
#' @return A data frame containing:
#' - Sim: the simulation run
#' - p0: the starting allele frequency after admixture
#' - pF: the final allele frequency of each simulation
#' - the parameter values used
#' @export

adsim_simulator <- function(
    ngens,              # Generations after admixture
    k,                  # Growth multiplier per generation
    l,                  # Generations before growth starts
    admix,              # Ancestry proportion from population 1
    m1,                 # Migration proportion from population 1
    m2,                 # Migration proportion from population 2
    ne,                 # Starting effective population size
    p1,                 # Allele frequency in population 1
    p2,                 # Allele frequency in population 2
    nsims,              # Number of simulations
    fixed_p0 = NULL,    # Optional fixed starting frequency
    max_n = 5000        # Maximum simulated population size
) {
  is_scalar_number <- function(x) {
    is.numeric(x) && length(x) == 1L && is.finite(x)
  }

  is_whole_number <- function(x) {
    is_scalar_number(x) && x == floor(x)
  }

  stopifnot(
    is_whole_number(ngens), ngens >= 0,
    is_whole_number(l), l >= 0, l <= ngens,
    is_whole_number(ne), ne > 0,
    is_whole_number(nsims), nsims > 0,
    is_whole_number(max_n), max_n > 0,
    ne <= max_n,
    is_scalar_number(k), k > 0,
    is_scalar_number(admix), admix >= 0, admix <= 1,
    is_scalar_number(m1), m1 >= 0, m1 <= 1,
    is_scalar_number(m2), m2 >= 0, m2 <= 1,
    m1 + m2 <= 1,
    is_scalar_number(p1), p1 >= 0, p1 <= 1,
    is_scalar_number(p2), p2 >= 0, p2 <= 1
  )

  if (!is.null(fixed_p0)) {
    stopifnot(
      is_scalar_number(fixed_p0),
      fixed_p0 >= 0,
      fixed_p0 <= 1
    )
  }

  # Set the starting frequency for each simulation.
  if (!is.null(fixed_p0)) {
    p0 <- rep(fixed_p0, nsims)
  } else {
    # Form the admixed population from the two source populations.
    raw_n1 <- admix * ne
    n1 <- floor(raw_n1) +
      rbinom(nsims, size = 1, prob = raw_n1 - floor(raw_n1))
    n2 <- ne - n1

    allele_count <- rbinom(nsims, size = 2 * n1, prob = p1) +
      rbinom(nsims, size = 2 * n2, prob = p2)

    p0 <- allele_count / (2 * ne)
  }

  current_p <- p0
  current_n <- ne

  # Simulate exactly ngens generations after the starting point.
  for (generation in seq_len(ngens)) {
    new_n <- if (generation > l) {
      floor(min(current_n * k, max_n))
    } else {
      ne
    }

    if (new_n < 1) {
      stop("Population size reached zero. Check k and ne.")
    }

    # Mix in migration; this has no effect when m1 and m2 are zero.
    mixed_p <- (1 - m1 - m2) * current_p +
      m1 * p1 + m2 * p2

    # Genetic drift: sample allele copies in the next generation.
    allele_count <- rbinom(nsims, size = 2 * new_n, prob = mixed_p)
    current_p <- allele_count / (2 * new_n)
    current_n <- new_n
  }

  data.frame(
    Sim = seq_len(nsims),
    p0 = p0,
    pF = current_p,
    ngens = ngens,
    k = k,
    l = l,
    admix = admix,
    m1 = m1,
    m2 = m2,
    ne = ne,
    p1 = p1,
    p2 = p2
  )
}

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
#' @param max_n Ceiling on the simulated census size, applied once population
#'   growth begins. Defaults to `Inf`, so size follows `k` and `l` without an
#'   arbitrary limit. Set it to a finite whole number to cap growth, for
#'   instance at a carrying capacity.
#' @details Each generation is drawn as a single vectorised binomial sample
#'   across all `nsims` runs, so cost scales with `nsims * ngens` and is
#'   effectively independent of `ne`: the sampler used by [stats::rbinom()]
#'   takes constant time in `size`, so an effective population size of 10^7
#'   costs no more per generation than one of 10^3.
#'
#'   Runs that reach an absorbing state are dropped from the sampling set,
#'   since their frequency can no longer change. Loss is absorbing when no
#'   migrant copies of the allele arrive (`m1 * p1 + m2 * p2 == 0`) and
#'   fixation is absorbing when no migrant non-copies arrive; with migration
#'   from a polymorphic source neither holds and every run stays active. This
#'   is an exact saving rather than an approximation, and it matters most for
#'   long runs at small `ne`, where most runs fix or are lost early.
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
    max_n = Inf         # Ceiling on simulated population size
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
    is.numeric(max_n), length(max_n) == 1L, !is.na(max_n), max_n > 0,
    is.infinite(max_n) || max_n == floor(max_n),
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

  # Census size for each generation, computed once rather than inside the
  # sampling loop.
  n_traj <- numeric(ngens)
  current_n <- ne
  for (generation in seq_len(ngens)) {
    new_n <- if (generation > l) {
      floor(min(current_n * k, max_n))
    } else {
      ne
    }

    if (new_n < 1) {
      stop("Population size reached zero. Check k and ne.")
    }

    n_traj[generation] <- new_n
    current_n <- new_n
  }

  # Set the starting frequency for each simulation.
  if (!is.null(fixed_p0)) {
    p0 <- rep(fixed_p0, nsims)
  } else {
    # Form the admixed population from the two source populations.
    raw_n1 <- admix * ne
    n1 <- floor(raw_n1) +
      stats::rbinom(nsims, size = 1, prob = raw_n1 - floor(raw_n1))
    n2 <- ne - n1

    allele_count <- stats::rbinom(nsims, size = 2 * n1, prob = p1) +
      stats::rbinom(nsims, size = 2 * n2, prob = p2)

    p0 <- allele_count / (2 * ne)
  }

  # Migration constants, hoisted out of the loop. The frequency entering the
  # drift step is mix_slope * p + mix_intercept.
  mix_slope <- 1 - m1 - m2
  mix_intercept <- (m1 * p1) + (m2 * p2)

  # A run can only stop changing if the state it reaches is absorbing. Loss is
  # absorbing when no migrant copies arrive, and fixation is absorbing when no
  # migrant non-copies arrive. With migration from a polymorphic source
  # neither holds and every run stays active.
  loss_absorbs <- (mix_intercept == 0)
  fixation_absorbs <- ((mix_slope + mix_intercept) == 1)
  prunable <- loss_absorbs || fixation_absorbs

  # Absorbed runs are removed from the sampling set, since their frequency can
  # no longer change. Checking costs a pass over the active runs, so it is
  # amortised over several generations rather than done every generation.
  check_every <- 8L

  pF <- numeric(nsims)        # final frequencies, filled in as runs absorb
  active_p <- p0              # frequencies of the runs still changing
  active_idx <- seq_len(nsims)

  for (generation in seq_len(ngens)) {
    two_n <- 2 * n_traj[generation]

    allele_count <- stats::rbinom(
      length(active_p), size = two_n,
      prob = (mix_slope * active_p) + mix_intercept
    )
    active_p <- allele_count / two_n

    if (prunable && (generation %% check_every == 0L)) {
      count_range <- range(allele_count)
      absorbed_present <-
        (loss_absorbs && count_range[1L] == 0) ||
        (fixation_absorbs && count_range[2L] == two_n)

      if (absorbed_present) {
        keep <- if (loss_absorbs && fixation_absorbs) {
          allele_count > 0 & allele_count < two_n
        } else if (loss_absorbs) {
          allele_count > 0
        } else {
          allele_count < two_n
        }

        pF[active_idx[!keep]] <- active_p[!keep]
        active_p <- active_p[keep]
        active_idx <- active_idx[keep]

        if (!length(active_idx)) break
      }
    }
  }

  pF[active_idx] <- active_p

  data.frame(
    Sim = seq_len(nsims),
    p0 = p0,
    pF = pF,
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

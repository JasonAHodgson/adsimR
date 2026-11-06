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
                  m1, m2, ne, p1, p2, nsims,
                  fixed_p0 = NULL,   # Optional fixed starting frequency
                  max_n = 5000) {    # Maximum simulated population size
  # Make every combination of parameter values.
  param_combos <- expand.grid(
    ngens = ngens, k = k, l = l, admix = admix,
    m1 = m1, m2 = m2, ne = ne, p1 = p1, p2 = p2, nsims = nsims,
    KEEP.OUT.ATTRS = FALSE
  )

  # Skip impossible combinations instead of stopping everything.
  valid <- param_combos$l <= param_combos$ngens &
    param_combos$m1 + param_combos$m2 <= 1
  if (!any(valid)) {
    stop("No valid combinations: check l <= ngens and m1 + m2 <= 1.")
  }
  if (any(!valid)) {
    warning(sum(!valid), " combination(s) skipped: l > ngens or m1 + m2 > 1.")
  }
  param_combos <- param_combos[valid, , drop = FALSE]

  message(nrow(param_combos), " scenario(s), ",
          sum(param_combos$nsims), " simulations in total.")

  # Run the simulator once per scenario; it checks each value itself.
  results_list <- lapply(seq_len(nrow(param_combos)), function(i) {
    pc <- param_combos[i, ]
    res <- adsim_simulator(
      ngens = pc$ngens, k = pc$k, l = pc$l, admix = pc$admix,
      m1 = pc$m1, m2 = pc$m2, ne = pc$ne, p1 = pc$p1, p2 = pc$p2,
      nsims = pc$nsims, fixed_p0 = fixed_p0, max_n = max_n
    )
    cbind(Scenario = i, res)
  })

  results <- do.call(rbind, results_list)
  rownames(results) <- NULL
  results
}


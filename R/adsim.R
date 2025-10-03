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


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
#' @param nsims Number of simulations to perform, or `"auto"`, the default, to
#'   choose a number for each scenario by simulating until the distribution is
#'   resolved to `precision`. A fixed number is used exactly as given; the
#'   precision it achieved is reported either way, so a fixed choice can be
#'   checked after the fact.
#' @param precision Target Monte Carlo half-width, on the 0 to 1 scale shared
#'   by allele frequencies and proportions. The default of 0.01 means every
#'   quantile and absorption proportion reported for a scenario is pinned to
#'   about one percentage point.
#' @param min_nsims Number of runs in the first batch when `nsims = "auto"`.
#' @param max_nsims Ceiling on the runs used for one scenario when
#'   `nsims = "auto"`. Reaching it stops the scenario with `converged = FALSE`
#'   in the summary rather than running without limit.
#' @param conf Confidence level for the Monte Carlo intervals.
#' @param observed Observed allele frequency the simulations will be tested
#'   against, used when `nsims = "auto"` to decide how deep to sample. Supply
#'   it and each scenario is run until the tail at that value is resolved to
#'   `tail_rse`; leave it `NULL` and the number of runs is decided by the shape
#'   of the distribution alone. It is worth supplying, because how many runs a
#'   test needs depends on where the observation falls in the null and barely
#'   at all on the scenario: an observation at the ninth decile is settled by a
#'   couple of hundred runs, one at the 99.9th percentile needs tens of
#'   thousands. The same value then goes to [adsim_pvalue()].
#' @param tail_p Smallest tail probability worth paying to resolve, used only
#'   when `observed` is given. Depth grows as 1 / `tail_p`, so an observation
#'   far out in the tail would otherwise demand runs without limit; the target
#'   is held here instead. It caps the cost rather than promising an accuracy,
#'   and `p_resolution` in the summary reports what was actually achieved.
#' @param tail_rse Relative standard error wanted for the tail probability at
#'   `observed`. The default of 0.2 asks for it to within about a fifth of
#'   itself.
#' @param fixed_p0 Optional fixed starting allele frequency, used for every run
#'   in place of forming generation 0 from the two parent populations.
#' @param max_n Ceiling on the simulated census size, applied once population
#'   growth begins. Defaults to `Inf`, so size follows `k` and `l` without an
#'   arbitrary limit. Set it to a finite whole number to cap growth, for
#'   instance at a carrying capacity.
#' @param keep_sims Whether to store the individual runs. When `FALSE` each
#'   scenario is summarised and its runs discarded, so `nsims` can be raised
#'   without the memory cost of keeping every run. Summary statistics are
#'   unaffected; only `$sims` and [as.data.frame()] need the runs.
#' @return An [adsim_result]: the conditions of each scenario, the individual
#'   runs unless `keep_sims = FALSE`, and a summary for each scenario covering
#'   the mean, variance and quantiles of final allele frequency, the proportion
#'   of runs fixed and lost, heterozygosity, and the exact expectations from
#'   [adsim_moments()]. Use [as.data.frame()] for the flat one-row-per-run
#'   format.
#' @details Generation 0 is formed by sampling `admix * ne` individuals from
#'   the first parent population and the remainder from the second. Each
#'   subsequent generation is produced in two steps. Migration acts first on
#'   the gene pool, following the continent-island model, so that the frequency
#'   entering the drift step is
#'   \deqn{p^{*} = (1 - m_1 - m_2)\,p + m_1 p_1 + m_2 p_2,}
#'   and genetic drift then follows, with \eqn{2N} gene copies drawn
#'   binomially at frequency \eqn{p^{*}}. Migration therefore changes allele
#'   frequency but not census size, which is governed by `k`, `l` and `max_n`
#'   alone.
#' @export

adsim <- function(ngens, k, l, admix,
                  m1, m2, ne, p1, p2, nsims = "auto",
                  fixed_p0 = NULL,   # Optional fixed starting frequency
                  max_n = Inf,       # Ceiling on simulated population size
                  keep_sims = TRUE,  # Store the individual runs?
                  precision = 0.01,  # Target Monte Carlo half-width
                  min_nsims = 500,   # First batch when choosing nsims
                  max_nsims = 1e5,   # Ceiling when choosing nsims
                  conf = 0.95,       # Confidence level for that half-width
                  observed = NULL,   # Observed frequency the depth targets
                  tail_p = 0.001,    # Smallest tail probability worth paying for
                  tail_rse = 0.2) {  # Relative error wanted for that tail
  choose_nsims <- is.character(nsims) &&
    length(nsims) == 1L && identical(nsims, "auto")
  if (!choose_nsims) {
    stopifnot(is.numeric(nsims), all(nsims > 0), all(nsims == floor(nsims)))
  }
  stopifnot(
    is.numeric(precision), length(precision) == 1L,
    precision > 0, precision < 1,
    is.numeric(min_nsims), length(min_nsims) == 1L, min_nsims >= 2,
    is.numeric(max_nsims), length(max_nsims) == 1L, max_nsims >= min_nsims,
    is.numeric(conf), length(conf) == 1L, conf > 0, conf < 1,
    is.numeric(tail_p), length(tail_p) == 1L, tail_p > 0, tail_p < 1,
    is.numeric(tail_rse), length(tail_rse) == 1L, tail_rse > 0
  )
  if (!is.null(observed)) {
    stopifnot(is.numeric(observed), length(observed) == 1L, !is.na(observed),
              observed >= 0, observed <= 1)
  }

  # Make every combination of parameter values. When nsims is chosen
  # adaptively it is not crossed with the rest; each scenario gets its own.
  param_combos <- expand.grid(
    ngens = ngens, k = k, l = l, admix = admix,
    m1 = m1, m2 = m2, ne = ne, p1 = p1, p2 = p2,
    nsims = if (choose_nsims) NA_real_ else nsims,
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

  if (choose_nsims) {
    message(nrow(param_combos), " scenario(s); choosing the number of runs ",
            "for each to reach a Monte Carlo half-width of ", precision, ".")
  } else {
    message(nrow(param_combos), " scenario(s), ",
            sum(param_combos$nsims), " simulations in total.")
  }

  # Run the simulator once per scenario; it checks each value itself.
  conditions <- vector("list", nrow(param_combos))
  summaries <- vector("list", nrow(param_combos))
  runs <- if (keep_sims) vector("list", nrow(param_combos)) else NULL

  for (i in seq_len(nrow(param_combos))) {
    pc <- param_combos[i, ]

    sim_args <- list(
      ngens = pc$ngens, k = pc$k, l = pc$l, admix = pc$admix,
      m1 = pc$m1, m2 = pc$m2, ne = pc$ne, p1 = pc$p1, p2 = pc$p2,
      fixed_p0 = fixed_p0, max_n = max_n
    )

    if (choose_nsims) {
      adaptive <- simulate_to_precision(
        sim_args, precision = precision, min_nsims = min_nsims,
        max_nsims = max_nsims, conf = conf, observed = observed,
        tail_p = tail_p, tail_rse = tail_rse
      )
      res <- data.frame(Sim = seq_along(adaptive$pF),
                        p0 = adaptive$p0, pF = adaptive$pF)
      pc$nsims <- length(adaptive$pF)
      reached <- adaptive$precision_reached
    } else {
      res <- do.call(adsim_simulator, c(sim_args, list(nsims = pc$nsims)))
      reached <- run_precision(res$pF, precision = precision, conf = conf)
      adaptive <- list(
        tail_depth = tail_depth(res$pF, observed = observed, tail_p = tail_p,
                                tail_rse = tail_rse)
      )
      adaptive$tails_resolved <- pc$nsims >= adaptive$tail_depth$required
    }

    expected <- adsim_moments(
      ngens = pc$ngens, k = pc$k, l = pc$l, admix = pc$admix,
      m1 = pc$m1, m2 = pc$m2, ne = pc$ne, p1 = pc$p1, p2 = pc$p2,
      fixed_p0 = fixed_p0, max_n = max_n
    )

    conditions[[i]] <- data.frame(
      scenario = i,
      ngens = pc$ngens, k = pc$k, l = pc$l, admix = pc$admix,
      m1 = pc$m1, m2 = pc$m2, ne = pc$ne, p1 = pc$p1, p2 = pc$p2,
      nsims = pc$nsims,
      fixed_p0 = if (is.null(fixed_p0)) NA_real_ else fixed_p0,
      max_n = max_n,
      n_final = expected$n_final,
      n_max = expected$n_max,
      n_harmonic = expected$n_harmonic
    )

    summaries[[i]] <- cbind(
      scenario = i,
      summarise_runs(res$p0, res$pF, expected),
      data.frame(
        mc_precision = reached$half_width,
        mc_precision_props = reached$half_width_props,
        mc_precision_quantiles = reached$half_width_quantiles,
        quantiles_at_atom = reached$at_atom,
        p_resolution = 1 / (pc$nsims + 1),
        tail_target = adaptive$tail_depth$target,
        tails_resolved = adaptive$tails_resolved,
        converged = reached$half_width <= precision &&
          adaptive$tails_resolved
      )
    )

    if (keep_sims) {
      runs[[i]] <- data.frame(
        scenario = i, sim = res$Sim, p0 = res$p0, pF = res$pF
      )
    }
  }

  new_adsim_result(
    conditions = do.call(rbind, conditions),
    sims = if (keep_sims) do.call(rbind, runs) else NULL,
    summary = do.call(rbind, summaries),
    call = match.call()
  )
}


#' Quantile levels reported by [adsim()]
#' @keywords internal
#' @noRd
reported_probs <- function() {
  c(0.025, 0.05, 0.25, 0.5, 0.75, 0.95, 0.975)
}

#' Wilson confidence half-width for a proportion
#'
#' @param successes Number of successes.
#' @param n Number of trials.
#' @param conf Confidence level.
#' @return The half-width of the Wilson interval.
#' @keywords internal
#' @noRd
wilson_half_width <- function(successes, n, conf = 0.95) {
  z <- stats::qnorm(1 - ((1 - conf) / 2))
  phat <- successes / n
  (z * sqrt((phat * (1 - phat) / n) + (z^2 / (4 * n^2)))) / (1 + (z^2 / n))
}

#' How precisely has a set of runs pinned down the distribution?
#'
#' Measures the Monte Carlo uncertainty of the quantities [adsim()] reports,
#' on the 0 to 1 scale they share.
#'
#' @param pF Numeric vector of final allele frequencies.
#' @param precision Target half-width, used to decide when a quantile is
#'   pinned to an atom.
#' @param conf Confidence level for the intervals.
#' @return A list with the worst half-width over all monitored quantities
#'   (`half_width`), its two parts (`half_width_props`, `half_width_quantiles`)
#'   and the number of quantile levels that sit on an atom (`at_atom`).
#' @details Proportions fixed and lost get a Wilson interval. Each reported
#'   quantile gets a distribution-free interval from order statistics, which
#'   makes no assumption about the shape of the distribution.
#'
#'   Quantile levels that fall on an atom are excluded and counted instead.
#'   Final frequencies pile up at 0 and 1, and when the cumulative distribution
#'   jumps across a level `q` the quantile at that level is that atom. Its
#'   order statistic interval then spans the gap to the next atom and stays
#'   wide however many runs are added: with half the runs lost and half fixed,
#'   the interval for the median runs from 0 to 1 at any number of runs, and
#'   resolving it would mean distinguishing a proportion from one half, which
#'   needs replicates without limit and answers nothing. The atom masses are
#'   controlled by the proportion criterion instead.
#' @keywords internal
#' @noRd
run_precision <- function(pF, precision = 0.01, conf = 0.95) {
  n <- length(pF)
  if (n < 2L) {
    return(list(half_width = Inf, half_width_props = Inf,
                half_width_quantiles = Inf, at_atom = 0L))
  }

  sorted <- sort(pF)
  alpha <- 1 - conf
  probs <- reported_probs()

  half_width_props <- max(
    wilson_half_width(sum(pF >= 1), n, conf),
    wilson_half_width(sum(pF <= 0), n, conf)
  )

  half_width_quantiles <- 0
  at_atom <- 0L
  for (q in probs) {
    value <- sorted[max(1L, ceiling(q * n))]
    cumulative_at <- mean(pF <= value)
    cumulative_below <- mean(pF < value)

    lower <- max(stats::qbinom(alpha / 2, n, q), 1L)
    upper <- min(stats::qbinom(1 - (alpha / 2), n, q) + 1L, n)
    this_half_width <- (sorted[upper] - sorted[lower]) / 2

    # The q-quantile is the smallest value whose cumulative probability
    # reaches q, so the atom spans the level when F(below) < q <= F(at).
    # The bound on the right must not be strict: an even split between two
    # atoms puts F(at) exactly on q.
    jumps_across <- (cumulative_at - cumulative_below) > precision &&
      cumulative_below < q && cumulative_at >= q

    # Only a level that both sits on an atom and remains wide is hopeless.
    # A level buried well inside an atom has a zero width interval and is
    # perfectly well determined, so it is kept.
    if (jumps_across && this_half_width > precision) {
      at_atom <- at_atom + 1L
    } else {
      half_width_quantiles <- max(half_width_quantiles, this_half_width)
    }
  }

  list(
    half_width = max(half_width_props, half_width_quantiles),
    half_width_props = half_width_props,
    half_width_quantiles = half_width_quantiles,
    at_atom = at_atom
  )
}

#' How many runs are needed to resolve the tail at an observed frequency?
#'
#' @param pF Numeric vector of final allele frequencies.
#' @param observed Observed allele frequency, or NULL for no requirement.
#' @param tail_p Smallest tail probability worth paying to resolve.
#' @param tail_rse Relative standard error wanted for the tail probability.
#' @return A list with the number of runs required and the tail probability
#'   the requirement was set for, both zero and `NA` when `observed` is NULL.
#' @details Estimating a tail probability \eqn{\alpha} from `n` runs gives it a
#'   relative standard error of \eqn{\sqrt{(1 - \alpha) / (n\alpha)}}, so
#'   resolving it to `tail_rse` needs
#'   \eqn{n = (1 - \alpha) / (\alpha \cdot \mathrm{tail\_rse}^2)} runs. That
#'   cost depends on \eqn{\alpha} alone and not at all on the distribution the
#'   probability came from, which is why the depth has to be keyed to the
#'   observation: it is where the observed frequency falls in the null, not
#'   the shape of the null, that decides how many runs are worth running. An
#'   observation at the ninth decile of its own null is settled by a couple of
#'   hundred runs whatever the scenario, while one at the 99.9th percentile
#'   needs tens of thousands.
#'
#'   The tail used is the smaller of the two, matching the two-sided default of
#'   [adsim_pvalue()]. An atom needs no special handling here: an observation
#'   at 1 when half the runs are fixed simply has an estimated tail of about a
#'   half, which asks for almost no depth, correctly, because no smaller
#'   p-value exists in that direction.
#'
#'   `tail_p` caps the cost rather than promising an accuracy. Depth grows as
#'   \eqn{1 / \alpha}, so an observation far out in the tail could demand runs
#'   without limit; the target is held at `tail_p` instead, and the p-value is
#'   then reported with `p_resolution` so its precision is visible.
#' @keywords internal
#' @noRd
tail_depth <- function(pF, observed = NULL, tail_p = 0.001, tail_rse = 0.2) {
  if (is.null(observed)) {
    return(list(required = 0L, target = NA_real_))
  }

  n <- length(pF)
  count <- min(sum(pF >= observed), sum(pF <= observed))
  alpha <- (count + 1) / (n + 1)
  target <- max(alpha, tail_p)

  list(
    required = ceiling((1 - target) / (target * (tail_rse^2))),
    target = target
  )
}

#' Simulate until the distribution is resolved to a target precision
#'
#' @param params Named list of simulator arguments, without `nsims`.
#' @param precision Target Monte Carlo half-width.
#' @param min_nsims Number of runs in the first batch.
#' @param max_nsims Ceiling on the total number of runs.
#' @param conf Confidence level.
#' @return A list with `p0`, `pF` and the precision reached.
#' @details Runs are added in batches. Monte Carlo half-widths shrink as
#'   `1 / sqrt(n)`, so after each batch the number still needed is estimated
#'   directly rather than by repeated doubling; growth is capped at eight fold
#'   per batch so that a noisy early estimate cannot overshoot badly.
#' @keywords internal
#' @noRd
simulate_to_precision <- function(params, precision = 0.01, min_nsims = 500L,
                                  max_nsims = 100000L, conf = 0.95,
                                  observed = NULL, tail_p = 0.001,
                                  tail_rse = 0.2) {
  p0 <- numeric(0)
  pF <- numeric(0)
  batch_size <- min(min_nsims, max_nsims)

  repeat {
    batch <- do.call(adsim_simulator, c(params, list(nsims = batch_size)))
    p0 <- c(p0, batch$p0)
    pF <- c(pF, batch$pF)
    done <- length(pF)

    reached <- run_precision(pF, precision = precision, conf = conf)
    depth <- tail_depth(pF, observed = observed, tail_p = tail_p,
                        tail_rse = tail_rse)

    shape_resolved <- reached$half_width <= precision
    tails_resolved <- done >= depth$required
    if ((shape_resolved && tails_resolved) || done >= max_nsims) break

    # half-widths scale as 1 / sqrt(n), so aim straight at the required total
    needed_for_shape <- if (shape_resolved) {
      done
    } else {
      ceiling(done * ((reached$half_width / precision)^2))
    }
    required <- max(needed_for_shape, depth$required)

    # Grow by at least a tenth so that progress is always made, but no more
    # than that beyond what the estimate asks for, or a scenario that is only
    # just short of the target overshoots it badly.
    batch_size <- min(max(required - done, ceiling(done / 10), 100),
                      8L * done, max_nsims - done)
    batch_size <- max(batch_size, 1L)
  }

  list(p0 = p0, pF = pF, precision_reached = reached, tail_depth = depth,
       tails_resolved = length(pF) >= depth$required)
}

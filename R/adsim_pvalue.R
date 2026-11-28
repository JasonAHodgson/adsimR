#' Test an observed allele frequency against simulated drift
#'
#' Compares an observed allele frequency with the distribution produced by
#' [adsim()], by counting how many simulated runs are at least as extreme. The
#' simulated runs are the null distribution, so no assumption is made about its
#' shape.
#'
#' @param x An [adsim_result] from [adsim()], simulated with
#'   `keep_sims = TRUE`.
#' @param observed Observed allele frequency, between 0 and 1.
#' @param alternative Direction of the test. `"two.sided"`, the default, makes
#'   no claim about the direction of selection. `"greater"` asks whether the
#'   allele is at a higher frequency than drift and admixture produce, and
#'   `"less"` the reverse.
#' @param conf Confidence level for the interval on the p-value.
#' @return A data frame with one row per scenario:
#' - scenario: the scenario the test refers to
#' - observed: the observed frequency tested
#' - alternative: the direction tested
#' - nsims: runs the null was built from
#' - n_beyond: runs at least as extreme as `observed`, ties included
#' - p_value: the p-value
#' - p_lower, p_upper: confidence interval for it
#' - p_resolution: the smallest p-value these runs could have produced
#' - resolved: whether `p_value` is above `p_resolution`, so is not merely a
#'   statement that the runs ran out
#' @details The p-value is \eqn{(r + 1) / (n + 1)}, where \eqn{r} counts the
#'   runs at least as extreme as `observed` and \eqn{n} is the number of runs.
#'   The added one is not a fudge: the observed value is itself one draw under
#'   the null, and without it a p-value of exactly zero is reported, which
#'   asserts an impossibility. The smallest value obtainable is therefore
#'   \eqn{1 / (n + 1)}, reported as `p_resolution`.
#'
#'   Ties count as extreme. Final frequencies pile up at 0 and 1, and excluding
#'   runs that equal the observed value would make the test anti-conservative
#'   exactly where those atoms sit.
#'
#'   The two-sided p-value is twice the smaller tail, capped at 1. Doubling is
#'   used rather than a region of equal density because the null is routinely
#'   asymmetric and often has two modes, which leaves no useful ordering of how
#'   extreme a value is on the two sides.
#'
#'   Where the null itself reaches an absorbing boundary the test has little
#'   power, and this is correct rather than a shortcoming. If drift alone fixes
#'   the allele in a tenth of runs, the smallest upper-tail p-value available
#'   is about a tenth, and an observation of fixation cannot be surprising. A
#'   normal approximation returns a small p-value for the same observation only
#'   by treating two atoms as a bell curve.
#' @seealso [adsim()] for the simulations, and [adsim_moments()] for the exact
#'   mean and standard deviation, which remain valid but should not be used to
#'   build a normal interval for a tail.
#' @export
#' @examples
#' sims <- adsim(ngens = 33, k = 1.1, l = 5, admix = 0.5, m1 = 0.01, m2 = 0.01,
#'               ne = 200, p1 = 0.5, p2 = 0.9, nsims = 2000)
#' adsim_pvalue(sims, observed = 0.95)

adsim_pvalue <- function(x, observed,
                         alternative = c("two.sided", "greater", "less"),
                         conf = 0.95) {
  alternative <- match.arg(alternative)

  if (!inherits(x, "adsim_result")) {
    stop("x must be an adsim_result, as returned by adsim().")
  }
  if (is.null(x$sims)) {
    stop("The individual runs are needed for this test. ",
         "Re-run adsim() with keep_sims = TRUE.")
  }
  stopifnot(
    is.numeric(observed), length(observed) == 1L, !is.na(observed),
    observed >= 0, observed <= 1,
    is.numeric(conf), length(conf) == 1L, conf > 0, conf < 1
  )

  alpha <- 1 - conf

  rows <- lapply(x$conditions$scenario, function(id) {
    runs <- x$sims$pF[x$sims$scenario == id]
    n <- length(runs)

    upper <- sum(runs >= observed)
    lower <- sum(runs <= observed)

    count <- switch(alternative,
      greater = upper,
      less = lower,
      two.sided = min(upper, lower)
    )

    p_value <- (count + 1) / (n + 1)
    # Clopper-Pearson bounds for the tail probability behind that count
    p_lower <- if (count == 0) 0 else stats::qbeta(alpha / 2, count, n - count + 1)
    p_upper <- if (count == n) 1 else stats::qbeta(1 - (alpha / 2), count + 1, n - count)

    if (alternative == "two.sided") {
      p_value <- min(1, 2 * p_value)
      p_lower <- min(1, 2 * p_lower)
      p_upper <- min(1, 2 * p_upper)
    }

    resolution <- (if (alternative == "two.sided") 2 else 1) / (n + 1)

    data.frame(
      scenario = id,
      observed = observed,
      alternative = alternative,
      nsims = n,
      n_beyond = count,
      p_value = p_value,
      p_lower = p_lower,
      p_upper = p_upper,
      p_resolution = min(1, resolution),
      resolved = p_value > min(1, resolution)
    )
  })

  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}

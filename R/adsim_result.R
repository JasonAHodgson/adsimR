#' Construct an adsim_result object
#'
#' @param conditions Data frame of simulation conditions, one row per scenario.
#' @param sims Data frame of per-run results, or NULL when runs were not kept.
#' @param summary Data frame of summary statistics, one row per scenario.
#' @param call The matched call that produced the object.
#' @return An object of class `adsim_result`.
#' @keywords internal
#' @noRd

new_adsim_result <- function(conditions, sims, summary, call = NULL) {
  structure(
    list(conditions = conditions, sims = sims, summary = summary),
    class = "adsim_result",
    created = Sys.time(),
    adsimR_version = as.character(utils::packageVersion("adsimR")),
    seed = NULL,
    call = call
  )
}

#' Summarise the runs of a single scenario
#'
#' @param p0 Numeric vector of starting allele frequencies.
#' @param pF Numeric vector of final allele frequencies.
#' @param expected One row data frame from [adsim_moments()].
#' @return A one row data frame of summary statistics.
#' @keywords internal
#' @noRd

summarise_runs <- function(p0, pF, expected) {
  n <- length(pF)
  mean_pF <- mean(pF)
  var_pF <- stats::var(pF)
  mean_p0 <- mean(p0)
  quantiles <- stats::quantile(
    pF, c(0.025, 0.05, 0.25, 0.5, 0.75, 0.95, 0.975), names = FALSE
  )

  data.frame(
    nsims = n,
    mean_p0 = mean_p0,
    sd_p0 = stats::sd(p0),
    mean_pF = mean_pF,
    var_pF = var_pF,
    sd_pF = sqrt(var_pF),
    se_mean_pF = sqrt(var_pF / n),
    se_sd_pF = sqrt(var_pF) / sqrt(2 * (n - 1)),
    q025 = quantiles[1], q05 = quantiles[2], q25 = quantiles[3],
    median = quantiles[4],
    q75 = quantiles[5], q95 = quantiles[6], q975 = quantiles[7],
    prop_fixed = mean(pF >= 1),
    prop_lost = mean(pF <= 0),
    prop_segregating = mean(pF > 0 & pF < 1),
    het_start = mean(2 * p0 * (1 - p0)),
    het_end = mean(2 * pF * (1 - pF)),
    fst = if (mean_p0 > 0 && mean_p0 < 1) {
      var_pF / (mean_p0 * (1 - mean_p0))
    } else {
      NA_real_
    },
    mean_expected = expected$mean,
    sd_expected = expected$sd,
    sd_ratio = if (expected$sd > 0) sqrt(var_pF) / expected$sd else NA_real_
  )
}

#' Coerce simulation output to a flat data frame
#'
#' Functions that accept the flat one-row-per-run format use this so that they
#' also accept an [adsim_result].
#'
#' @param x An `adsim_result` or a data frame.
#' @return A data frame with one row per run.
#' @keywords internal
#' @noRd

as_run_frame <- function(x) {
  if (inherits(x, "adsim_result")) as.data.frame(x) else x
}

#' Simulation results for one or more sets of conditions
#'
#' `adsim()` returns an object of class `adsim_result`. It holds the conditions
#' of each scenario once, rather than repeating them on every simulated run, and
#' carries summary statistics computed for each scenario.
#'
#' @section Components:
#' \describe{
#'   \item{conditions}{One row per scenario: the simulation parameters, plus
#'     `n_final`, `n_max` and `n_harmonic`, the final, largest and harmonic mean
#'     census sizes. `n_harmonic` is the size that governs drift, because drift
#'     accumulates as a sum of \eqn{1 / 2N} terms, and it can be far below
#'     `n_final` in a growing population.}
#'   \item{sims}{One row per run, with columns `scenario`, `sim`, `p0` and `pF`.
#'     `NULL` when `adsim()` was called with `keep_sims = FALSE`.}
#'   \item{summary}{One row per scenario. See below.}
#' }
#'
#' @section Summary columns:
#' \describe{
#'   \item{mean_p0, sd_p0, mean_pF, var_pF, sd_pF}{Location and spread of the
#'     starting and final allele frequencies.}
#'   \item{se_mean_pF, se_sd_pF}{Monte Carlo error of the mean and of the
#'     standard deviation, for judging whether `nsims` is large enough.
#'     `se_sd_pF` assumes approximate normality and understates the error for
#'     the strongly skewed distributions that arise near fixation.}
#'   \item{q025 to q975}{Empirical quantiles of the final allele frequency.
#'     Prefer these to intervals built from `sd_pF`: the distribution is often
#'     skewed, and U-shaped once runs begin to fix or be lost.}
#'   \item{prop_fixed, prop_lost, prop_segregating}{Proportion of runs ending
#'     fixed, lost, or still segregating.}
#'   \item{het_start, het_end}{Mean expected heterozygosity \eqn{2p(1 - p)} at
#'     the founding generation and at the end.}
#'   \item{fst}{Standardised variance among runs,
#'     `var_pF / (mean_p0 * (1 - mean_p0))`. It runs from 0 to 1 as runs
#'     differentiate, reaching 1 once every run has fixed or been lost. It
#'     equals `1 - het_end / het_start` exactly when the mean is unchanged (no
#'     migration) and the starting frequency is the same in every run, which
#'     needs `fixed_p0`: a sampled founding generation leaves variance in `p0`
#'     and the identity becomes approximate. Even then the two sides differ by
#'     the Bessel correction, since `var_pF` divides by `nsims - 1` while
#'     `het_end` is a plain mean over runs.
#'     Under migration the identity does not hold and heterozygosity can
#'     increase, because migration pulls the frequency toward an intermediate
#'     value; that is why the two heterozygosities are reported raw rather than
#'     as a proportion retained.}
#'   \item{mean_expected, sd_expected, sd_ratio}{The exact mean and standard
#'     deviation from [adsim_moments()], and the ratio of simulated to expected
#'     standard deviation. A ratio near 1 confirms the simulation behaved as
#'     the conditions imply. It is a diagnostic, not a test, and is noisy when
#'     `nsims` is small.}
#' }
#'
#' @name adsim_result
NULL

#' @describeIn adsim_result Print a compact overview.
#' @param x An `adsim_result`.
#' @param ... Ignored.
#' @return `print()` returns `x` invisibly.
#' @export

print.adsim_result <- function(x, ...) {
  cat("<adsim_result> ", nrow(x$conditions), " scenario(s), ",
      format(sum(x$conditions$nsims), big.mark = ","), " simulated run(s)\n",
      sep = "")
  cat("  created ", format(attr(x, "created"), "%Y-%m-%d %H:%M:%S"),
      " with adsimR ", attr(x, "adsimR_version"), "\n", sep = "")

  varying <- names(x$conditions)[
    vapply(x$conditions, function(column) length(unique(column)) > 1L,
           logical(1))
  ]
  varying <- setdiff(varying, c("scenario", "n_final", "n_max", "n_harmonic"))
  cat("  varying: ",
      if (length(varying)) paste(varying, collapse = ", ") else "nothing",
      "\n", sep = "")
  cat("  runs stored: ", if (is.null(x$sims)) "no" else "yes", "\n", sep = "")

  shown <- c("scenario", "mean_pF", "sd_pF", "prop_fixed", "prop_lost")
  print(utils::head(x$summary[, shown, drop = FALSE], 10L))
  if (nrow(x$summary) > 10L) {
    cat("  ... ", nrow(x$summary) - 10L, " further scenario(s)\n", sep = "")
  }
  invisible(x)
}

#' @describeIn adsim_result Return the per-scenario summary table.
#' @param object An `adsim_result`.
#' @export

summary.adsim_result <- function(object, ...) {
  object$summary
}

#' @describeIn adsim_result Expand to the flat one-row-per-run format, joining
#'   each run to the conditions that produced it.
#' @param row.names Passed to the data frame method; ignored here.
#' @param optional Passed to the data frame method; ignored here.
#' @export

as.data.frame.adsim_result <- function(x, row.names = NULL, optional = FALSE,
                                       ...) {
  if (is.null(x$sims)) {
    stop("Runs were not stored. Re-run adsim() with keep_sims = TRUE.")
  }
  flat <- merge(x$sims, x$conditions, by = "scenario", sort = FALSE)
  flat <- flat[order(flat$scenario, flat$sim), , drop = FALSE]
  rownames(flat) <- NULL
  flat
}

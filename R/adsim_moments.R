#' Census size in each simulated generation
#'
#' Reproduces the deterministic census trajectory used by [adsim_simulator()].
#'
#' @inheritParams adsim_simulator
#' @return A numeric vector of length `ngens` giving the census size in each
#'   generation after the founding generation.
#' @keywords internal
#' @noRd

census_trajectory <- function(ngens, k, l, ne, max_n = Inf) {
  trajectory <- numeric(ngens)
  current_n <- ne
  for (generation in seq_len(ngens)) {
    current_n <- if (generation > l) {
      floor(min(current_n * k, max_n))
    } else {
      ne
    }
    trajectory[generation] <- current_n
  }
  trajectory
}

#' Expected mean and standard deviation of final allele frequency
#'
#' Computes the mean and standard deviation of the final allele frequency for a
#' set of simulation conditions exactly, without simulating. Use it to predict
#' the outcome of a scenario, to check a simulation, or to avoid simulation
#' entirely when only the first two moments are needed.
#'
#' @inheritParams adsim_simulator
#' @return A one row data frame with columns:
#' - mean: expected final allele frequency
#' - sd: standard deviation of final allele frequency across runs
#' - var: variance of final allele frequency across runs
#' - n_final: census size in the final generation
#' - n_max: largest census size reached
#' - n_harmonic: harmonic mean census size across generations
#' @details Genetic drift is binomial sampling, whose variance depends only on
#'   the frequency entering the drift step. That frequency is a linear function
#'   of the current frequency, so the recursion over the first two moments
#'   closes and no distributional assumption is needed. Writing
#'   \eqn{a = 1 - m_1 - m_2} and \eqn{b = m_1 p_1 + m_2 p_2}, so that the
#'   migration step is \eqn{p^{*} = ap + b}, one generation at census size
#'   \eqn{N} takes \eqn{(\mu, \sigma^2)} to
#'   \deqn{\mu' = a\mu + b, \qquad
#'         \sigma^{2\prime} = a^{2}\sigma^{2}\left(1 - \frac{1}{2N}\right) +
#'           \frac{\mu'(1 - \mu')}{2N}.}
#'   This is iterated over the census trajectory implied by `k`, `l` and
#'   `max_n`. The founding generation contributes
#'   \eqn{E[p_0] = a_{dm} p_1 + (1 - a_{dm}) p_2} and
#'   \deqn{\mathrm{Var}[p_0] = \frac{A p_1 (1 - p_1) + (n_e - A) p_2 (1 - p_2)}
#'                                  {2 n_e^{2}}
#'         + \frac{(p_1 - p_2)^{2} f (1 - f)}{n_e^{2}},}
#'   where \eqn{A = a_{dm} n_e} and \eqn{f} is its fractional part; the second
#'   term is the extra variance from rounding the founder count stochastically.
#'   When `fixed_p0` is supplied the starting variance is zero.
#'
#'   With no migration and constant census size this reduces to the familiar
#'   drift decay \eqn{p_0 (1 - p_0) [1 - (1 - 1 / 2N)^{t}]}.
#'
#'   Only the mean and variance are exact; the shape of the distribution is not
#'   given. Final frequencies are often strongly skewed, and U-shaped once an
#'   appreciable proportion of runs have fixed or been lost, so a normal
#'   interval built from this standard deviation will misstate tail
#'   probabilities. Use simulated quantiles for those.
#' @seealso [adsim()], which reports these alongside the simulated values.
#' @export
#' @examples
#' adsim_moments(ngens = 33, k = 1.1, l = 5, admix = 0.5,
#'               m1 = 0.01, m2 = 0.01, ne = 200, p1 = 0.5, p2 = 0.9)

adsim_moments <- function(ngens, k, l, admix, m1, m2, ne, p1, p2,
                          fixed_p0 = NULL, max_n = Inf) {
  stopifnot(
    is.numeric(ngens), length(ngens) == 1L, ngens >= 0, ngens == floor(ngens),
    is.numeric(l), length(l) == 1L, l >= 0, l <= ngens,
    is.numeric(ne), length(ne) == 1L, ne > 0, ne == floor(ne),
    is.numeric(k), length(k) == 1L, k > 0,
    is.numeric(admix), length(admix) == 1L, admix >= 0, admix <= 1,
    is.numeric(m1), length(m1) == 1L, m1 >= 0, m1 <= 1,
    is.numeric(m2), length(m2) == 1L, m2 >= 0, m2 <= 1,
    m1 + m2 <= 1,
    is.numeric(p1), length(p1) == 1L, p1 >= 0, p1 <= 1,
    is.numeric(p2), length(p2) == 1L, p2 >= 0, p2 <= 1,
    is.numeric(max_n), length(max_n) == 1L, max_n > 0, ne <= max_n
  )
  if (!is.null(fixed_p0)) {
    stopifnot(is.numeric(fixed_p0), length(fixed_p0) == 1L,
              fixed_p0 >= 0, fixed_p0 <= 1)
  }

  trajectory <- census_trajectory(ngens, k, l, ne, max_n)

  if (is.null(fixed_p0)) {
    founders <- admix * ne
    fraction <- founders - floor(founders)
    mu <- (admix * p1) + ((1 - admix) * p2)
    sigma2 <-
      ((founders * p1 * (1 - p1)) + ((ne - founders) * p2 * (1 - p2))) /
        (2 * ne^2) +
      (((p1 - p2)^2) * fraction * (1 - fraction)) / (ne^2)
  } else {
    mu <- fixed_p0
    sigma2 <- 0
  }

  mix_slope <- 1 - m1 - m2
  mix_intercept <- (m1 * p1) + (m2 * p2)

  for (generation in seq_len(ngens)) {
    two_n <- 2 * trajectory[generation]
    mu <- (mix_slope * mu) + mix_intercept
    sigma2 <- ((mix_slope^2) * sigma2 * (1 - (1 / two_n))) +
      ((mu * (1 - mu)) / two_n)
  }

  data.frame(
    mean = mu,
    sd = sqrt(sigma2),
    var = sigma2,
    n_final = if (ngens > 0) trajectory[ngens] else ne,
    n_max = if (ngens > 0) max(trajectory) else ne,
    n_harmonic = if (ngens > 0) ngens / sum(1 / trajectory) else ne
  )
}

#' Compute selection coefficients from allele frequency scenarios
#'
#' This function takes a data frame or a path to an input file containing allele
#' frequency scenarios and computes the selection coefficients based on the
#' provided data.
#'
#' @param infile Either a data frame or path to input file, containing the
#' following columns:
#' - g generations
#' - q: initial allele frequency
#' - qx: final allele frequency
#' @param outfile Path to save output table. This will include added column s
#'   giving the selection coefficient.
#' @param recessive Logical, if TRUE assumes recessive model, if FALSE assumes
#'   dominant model. Default is FALSE.
#' @details Selection is modelled as viability selection at a single diallelic
#'   locus under Hardy-Weinberg proportions. Writing the focal allele frequency
#'   as \eqn{q}, with \eqn{p = 1 - q}, the relative genotype fitnesses are
#'   \eqn{w_{pp} = 1}, \eqn{w_{pq} = 1 + hs} and \eqn{w_{qq} = 1 + s}, where
#'   \eqn{h} is the dominance coefficient of the focal allele. Setting
#'   `recessive = TRUE` gives \eqn{h = 0} (the focal allele is fully recessive,
#'   so selection acts only on \eqn{qq} homozygotes) and `recessive = FALSE`
#'   gives \eqn{h = 1} (fully dominant, so heterozygotes have the same fitness
#'   as \eqn{qq} homozygotes). The allele frequency after one generation is
#'   \deqn{q' = \frac{q\,[\,p(1 + hs) + q(1 + s)\,]}{1 + 2hspq + sq^{2}}.}
#'   For a given observed frequency change a dominant allele therefore requires
#'   a smaller selection coefficient than a recessive one, because selection
#'   also acts on heterozygotes.
#'
#'   The coefficient is found by grid search over \eqn{s} in steps of
#'   \eqn{10^{-4}}, returning the smallest value whose projected frequency
#'   matches `qx` to two decimal places. The recessive case reproduces the
#'   original `s_co.pl` Perl implementation in `data-raw/perl_base_code`.
#' @return A data frame with selection coefficients, containing the columns:
#' - g: generations
#' - q: initial allele frequency
#' - qx: final allele frequency
#' - s: selection coefficient
#' @export
#' @examples
#' result <- compute_selection_coefficients(infile = system.file("perl_base_code/Selection_scenarios.txt", package = "adsimR"),
#' outfile = paste0(tempdir(), "/output.txt"))
compute_selection_coefficients <- function(infile, outfile = NULL,
                                           recessive = FALSE,
                                           s_max = 1) {
  # Read input: a data frame or a file path.
  if (is.data.frame(infile)) {
    dat <- infile
  } else if (is.character(infile) && length(infile) == 1L) {
    if (!file.exists(infile)) stop("ERROR: Cannot find ", infile)
    dat <- read.table(infile, header = TRUE)
  } else {
    stop("ERROR: 'infile' must be a data frame or a single file path")
  }

  if (!all(c("g", "q", "qx") %in% colnames(dat))) {
    stop("Input must contain columns: g, q, qx")
  }
  stopifnot(
    is.numeric(dat$g), all(dat$g >= 0), all(dat$g == floor(dat$g)),
    is.numeric(dat$q), all(dat$q >= 0 & dat$q <= 1),
    is.numeric(dat$qx), all(dat$qx >= 0 & dat$qx <= 1),
    is.numeric(s_max), length(s_max) == 1L, s_max > 0
  )

  # One generation of selection on the allele (frequency q).
  a_prime <- function(q, s, recessive) {
    p <- 1 - q
    if (recessive) {
      # Only qq homozygotes have fitness 1 + s.
      (p * q + q^2 * (1 + s)) / (1 + s * q^2)
    } else {
      # Carriers (pq and qq) have fitness 1 + s.
      (p * q * (1 + s) + q^2 * (1 + s)) / (1 + s * (q^2 + 2 * p * q))
    }
  }

  # Repeat for g generations.
  get_qx <- function(q, g, s, recessive) {
    for (i in seq_len(g)) q <- a_prime(q, s, recessive)
    q
  }

  # Find the s where the predicted final frequency equals the observed qx.
  get_s <- function(q, qx, g, recessive) {
    f <- function(s) get_qx(q, g, s, recessive) - qx
    lo <- f(0)
    hi <- f(s_max)
    if (lo == 0) return(0)
    if (lo > 0 || hi < 0) return(NA_real_)  # qx not reachable with 0 <= s <= s_max
    uniroot(f, c(0, s_max), tol = 1e-8)$root
  }

  dat$s <- mapply(get_s, dat$q, dat$qx, dat$g,
                  MoreArgs = list(recessive = recessive))

  if (!is.null(outfile)) {
    write.table(dat, outfile, quote = FALSE, sep = "\t", row.names = FALSE)
  }
  dat
}


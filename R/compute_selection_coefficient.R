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
#' result <-
#'       compute_selection_coefficients(infile = system.file("extdata",
#'                                                           "Selection_scenarios.txt",
#'                                                            package = "adsimR"),
#'                                    outfile = paste0(tempdir(), "/output.txt"))

compute_selection_coefficients <- function(infile, outfile, recessive = FALSE) {
  if(is.data.frame(infile)) {
    dat <- infile
  } else if (is.character(infile)) {
    if(!file.exists(infile)) stop(paste("ERROR: Cannot find", infile))
    dat <- utils::read.table(infile, header = TRUE)
  } else {
    stop("ERROR: 'infile' must be either a data frame or a file path (character string)")
  }

  if (!all(c("g", "q", "qx") %in% colnames(dat))) {
    stop("Input file must contain columns: g, q, qx")
  }

  # check that q and qx are rounded to two decimal places
  if (any(round(dat$q, 2) != dat$q) || any(round(dat$qx, 2) != dat$qx)) {
    message("Rounding q and qx to two decimal places for consistency.")
    dat$q <- round(dat$q, 2)
    dat$qx <- round(dat$qx, 2)
  }

  # One generation of viability selection at a single diallelic locus.
  # Relative genotype fitnesses, with `allele` the focal allele q:
  #   w(pp) = 1, w(pq) = 1 + h * s, w(qq) = 1 + s
  # where h is the dominance coefficient of the focal allele: h = 0 makes the
  # focal allele fully recessive, h = 1 makes it fully dominant. The mean
  # fitness must weight heterozygotes by their Hardy-Weinberg frequency 2pq.
  a_prime <- function(allele, s, recessive) {
    h <- if (recessive) 0 else 1
    p <- 1 - allele
    w_bar <- 1 + (2 * h * s * p * allele) + (s * allele^2)
    (allele * ((p * (1 + h * s)) + (allele * (1 + s)))) / w_bar
  }

  get_qx <- function(q, g, s, recessive) {
    for (i in seq_len(g)) {
      q <- a_prime(q, s, recessive)
    }
    q
  }

  get_s <- function(q, qx, g, recessive) {
    s_vals <- seq(-0.999, 0.9999, by = 0.0001)
    for (s in s_vals) {
      testq <- get_qx(q, g, s, recessive)
      if (round(testq, 2) == qx) return(s)
    }
    message("Selection coefficient is greater than 0.9999 or less ",
            "than -0.9999, check your input values.")
    return(NA)
  }

  dat$s <- mapply(get_s, dat$q, dat$qx, dat$g, recessive)

  utils::write.table(dat, outfile, quote = FALSE, sep = "\t", row.names = FALSE)

  return(dat)
}

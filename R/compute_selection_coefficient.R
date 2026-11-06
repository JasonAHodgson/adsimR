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

  a_prime <- function(allele, s, recessive) {
    if(recessive){
      p <- 1 - allele
      ((p * allele) + allele^2 + s * allele^2) / (1 + s * allele^2) # increasing denominator in proportion
    } else {
      q <- 1 - allele
      ((allele*q) + s*(allele*q) + allele^2 + s*(allele^2)) / (1 + s * allele^2 + s * (allele * q))
    }
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

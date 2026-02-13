#' Compute selection coefficients from allele frequency scenarios
#'
#' @param infile Path to input file (columns: g, q, qx).
#' @param outfile Path to save output table (with added column s).
#' @param recessive Logical, if TRUE assumes recessive model, if FALSE assumes
#'   dominant model. Default is FALSE.
#' @return Data frame with selection coefficients.
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

  a_prime <- function(allele, s, recessive) {
    if(recessive){
      p <- 1 - allele
      ((p * allele) + allele^2 + s * allele^2) / (1 + s * allele^2) # increasing denominator in proportion
    } else {
      q <- 1 - allele
      ((allele*q) + s*(allele*q) + allele^2 + s*(allele^2) / (1 + s * allele^2 + s * (allele * q)))
    }
  }

  get_qx <- function(q, g, s, recessive) {
    for (i in seq_len(g)) {
      q <- a_prime(q, s, recessive)
    }
    q
  }

  get_s <- function(q, qx, g, recessive) {
    s_vals <- seq(0.0001, 0.9999, by = 0.0001)
    for (s in s_vals) {
      testq <- get_qx(q, g, s, recessive)
      if (round(testq, 2) == qx) return(s)
    }
    return(NA)
  }

  dat$s <- mapply(get_s, dat$q, dat$qx, dat$g, recessive)

  utils::write.table(dat, outfile, quote = FALSE, sep = "\t", row.names = FALSE)

  return(dat)
}

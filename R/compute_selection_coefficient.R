#' Compute selection coefficients from allele frequency scenarios
#'
#' @param infile Path to input file (tab-delimited, with columns: g, q, qx).
#' @param outfile Path to save output table (with added column s).
#' @return Data frame with selection coefficients.
#' @examples
#' result <- compute_selection_coefficients(infile = system.file("perl_base_code/Selection_scenarios.txt", package = "adsimR"),
#' outfile = paste0(tempdir(), "/output.txt"))

compute_selection_coefficients <- function(infile, outfile) {
  if (!file.exists(infile)) stop(paste("ERROR: Cannot find", infile))

  dat <- read.table(infile, header = TRUE)

  if (!all(c("g", "q", "qx") %in% colnames(dat))) {
    stop("Input file must contain columns: g, q, qx")
  }

  q_prime <- function(q, s) {
    p <- 1 - q
    ((p * q) + q^2 + s * q^2) / (1 + s * q^2)
  }

  get_qx <- function(q, g, s) {
    for (i in seq_len(g)) {
      q <- q_prime(q, s)
    }
    q
  }

  get_s <- function(q, qx, g) {
    s_vals <- seq(0.0001, 0.9999, by = 0.0001)
    for (s in s_vals) {
      testq <- get_qx(q, g, s)
      if (round(testq, 2) == qx) return(s)
    }
    return(NA)
  }

  dat$s <- mapply(get_s, dat$q, dat$qx, dat$g)

  write.table(dat, outfile, quote = FALSE, sep = "\t", row.names = FALSE)

  return(dat)
}

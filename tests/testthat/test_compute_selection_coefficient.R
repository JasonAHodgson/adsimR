test_that("test valid input", {
  # create basic dataframe
  df <-
    data.frame(c = c(1, 2, 3), b = c(0.1, 0.2, 0.3), a = c(0.15, 0.25, 0.35))
  # save as text file
  input_file <- paste0(tempfile("test_input.txt"))
  write.table(df, input_file, sep = "\t", row.names = FALSE, quote = FALSE)
  expect_error(compute_selection_coefficients(
    infile = input_file, recessive = TRUE,
    outfile = paste0(tempdir(), "/output.txt")
  ), "Input file must contain columns: g, q, qx")
})

test_that("test valid csv input", {
  # create basic dataframe
  df <-
    data.frame(g = c(1, 2, 3), q = c(0.1, 0.2, 0.3), qx = c(0.15, 0.25, 0.35))
  # create a csv path
  input_file <- paste0(tempfile("test_input.csv"))
  write.table(df,
    input_file,
    col.names = c("g", "q", "qx"),
    row.names = FALSE,
    quote = FALSE
  )
  results <-
    compute_selection_coefficients(
      infile = input_file,
      recessive = TRUE,
      outfile = paste0(tempdir(), "/output.txt")
    )
  expect_true("s" %in% colnames(results))
  expect_equal(nrow(results), 3)
})

test_that("results are compatible with perl script", {
  path <- system.file("extdata", "Selection_scenarios.txt", package = "adsimR")
  results_adsimR <-
    compute_selection_coefficients(
      infile = path,
      recessive = TRUE,
      outfile = paste0(tempdir(), "/output.txt")
    )
  path2 <- system.file("extdata",
    "Selection_scenarios_coefficients.txt",
    package = "adsimR"
  )
  results_perl <- read.table(path2, header = TRUE)
  expect_equal(results_adsimR, results_perl)
})

test_that("the dominance effect reverses with starting allele frequency", {
  # While a beneficial allele is rare it spreads faster when dominant, because
  # selection acts on heterozygotes as well as on homozygotes, so a smaller s
  # accounts for a given change. As the allele approaches fixation the opposing
  # alleles are sheltered in heterozygotes of equal fitness and progress
  # stalls, so there a recessive allele advances faster and needs the smaller
  # s. The crossover lies near q = 0.5.
  s_of <- function(dat, recessive) {
    compute_selection_coefficients(
      infile = dat,
      outfile = file.path(tempdir(), "dominance.txt"),
      recessive = recessive
    )$s
  }

  rare <- data.frame(g = c(30, 40), q = c(0.05, 0.05), qx = c(0.20, 0.30))
  expect_true(all(s_of(rare, TRUE) > s_of(rare, FALSE)))

  common <- data.frame(g = c(33, 43), q = c(0.67, 0.67), qx = c(0.92, 0.92))
  expect_true(all(s_of(common, FALSE) > s_of(common, TRUE)))
})

test_that("dominant model recovers a known selection coefficient", {
  # Independent implementation of the one-locus diploid recursion, used here as
  # ground truth: w(pp) = 1, w(pq) = 1 + h * s, w(qq) = 1 + s.
  q_next <- function(q, s, h) {
    p <- 1 - q
    (q * ((p * (1 + h * s)) + (q * (1 + s)))) /
      (1 + (2 * h * s * p * q) + (s * q^2))
  }
  project <- function(q, s, g, h) {
    for (i in seq_len(g)) q <- q_next(q, s, h)
    q
  }

  s_true <- 0.05
  scenarios <- data.frame(g = c(20, 30, 40), q = c(0.20, 0.35, 0.50))
  scenarios$qx <-
    round(mapply(project, scenarios$q, s_true, scenarios$g, 1), 2)

  results <- compute_selection_coefficients(
    infile = scenarios,
    outfile = file.path(tempdir(), "dominant.txt"),
    recessive = FALSE
  )

  # The grid search returns the smallest s whose projection matches qx to two
  # decimal places, so s is recovered to the width of that rounding window.
  expect_true(all(abs(results$s - s_true) < 0.01))
  # The returned s must reproduce qx under the analytic recursion.
  expect_equal(
    round(mapply(project, results$q, results$s, results$g, 1), 2),
    results$qx
  )
})

test_that("recessive model matches the analytic recursion at h = 0", {
  q_next <- function(q, s, h) {
    p <- 1 - q
    (q * ((p * (1 + h * s)) + (q * (1 + s)))) /
      (1 + (2 * h * s * p * q) + (s * q^2))
  }
  project <- function(q, s, g, h) {
    for (i in seq_len(g)) q <- q_next(q, s, h)
    q
  }

  s_true <- 0.08
  scenarios <- data.frame(g = c(25, 50), q = c(0.30, 0.45))
  scenarios$qx <-
    round(mapply(project, scenarios$q, s_true, scenarios$g, 0), 2)

  results <- compute_selection_coefficients(
    infile = scenarios,
    outfile = file.path(tempdir(), "recessive.txt"),
    recessive = TRUE
  )

  expect_true(all(abs(results$s - s_true) < 0.01))
})

test_that("rounding q and qx", {
  # create basic dataframe
  df <-
    data.frame(g = c(1, 2, 3),
               q = c(0.1, 0.2, 0.3),
               qx = c(0.1519, 0.2545, 0.3556))
  expect_message(
    compute_selection_coefficients(
      infile = df,
      recessive = TRUE,
      outfile = paste0(
        tempdir(),
        "/output.txt"
      )
    ),
    "Rounding q and qx to two decimal places for consistency."
  )
})

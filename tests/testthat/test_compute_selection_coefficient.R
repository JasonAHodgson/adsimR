test_that("test valid input", {
  #create basic dataframe
  df <- data.frame(c = c(1, 2, 3), b = c(0.1, 0.2, 0.3), a = c(0.15, 0.25, 0.35))
  #save as text file
  input_file <- paste0(tempfile("test_input.txt"))
  write.table(df, input_file, sep = "\t", row.names = FALSE, quote = FALSE)
  expect_error(compute_selection_coefficients(infile = input_file, recessive = TRUE,
  outfile = paste0(tempdir(), "/output.txt")), "Input file must contain columns: g, q, qx")
})

test_that("test valid csv input", {
  #create basic dataframe
  df <- data.frame(g = c(1, 2, 3), q = c(0.1, 0.2, 0.3), qx = c(0.15, 0.25, 0.35))
  #create a csv path
  input_file <- paste0(tempfile("test_input.csv"))
  write.table(df, input_file, col.names = c("g", "q", "qx"), row.names = FALSE, quote = FALSE)
  results <- compute_selection_coefficients(infile = input_file, recessive = TRUE, outfile = paste0(tempdir(), "/output.txt"))
  expect_true("s" %in% colnames(results))
  expect_equal(nrow(results), 3)
})

test_that("results are compatible with perl script", {
  path <- system.file("/perl_base_code/Selection_scenarios.txt", package = "adsimR")
  results_adsimR <- compute_selection_coefficients(infile = path, recessive = TRUE, outfile = paste0(tempdir(), "/output.txt"))
  path2 <- system.file("/perl_base_code/Selection_scenarios_coefficients.txt", package = "adsimR")
  results_perl <- read.table(path2, header = TRUE)
  expect_equal(results_adsimR, results_perl)
})

test_that("dominant is less than recesive",{
  path <- system.file("/perl_base_code/Selection_scenarios.txt", package = "adsimR")
  results_recessive <- compute_selection_coefficients(infile = path, outfile = paste0(tempdir(), "/output.txt"), recessive = TRUE)
  results_dominant <- compute_selection_coefficients(infile = path, outfile = paste0(tempdir(), "/output.txt"), recessive = FALSE)
  expect_equal(results_recessive$s > results_dominant$s, rep(TRUE, 15))
})

test_that("adsim returns the right number of rows and >= 12 columns", {
  # parameters
  ngens <- c(20, 30)
  k     <- c(1, 1.05)
  l     <- 1
  admix <- c(0.5, 0.6)
  m1    <- 0
  m2    <- 0
  ne    <- c(100, 500)
  p1    <- 1
  p2    <- 0
  nsims <- 1000

  out <- adsim(ngens, k, l, admix, m1, m2, ne, p1, p2, nsims)

  # expected rows = sum(nsims across all parameter combinations (works for scalar or vector nsims)
  grid <- expand.grid(
    ngens = ngens, k = k, l = l, admix = admix,
    m1 = m1, m2 = m2, ne = ne, p1 = p1, p2 = p2, nsims = nsims
  )
  expected_rows <- sum(grid$nsims)

  expect_equal(nrow(out), expected_rows)
  expect_true(ncol(out) >= 12L)
  expect_true(all(c("Sim","p0","pF",
                    "ngens","k","l","admix","m1","m2","ne","p1","p2") %in% names(out)))
})

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
  path <- system.file("extdata","Selection_scenarios.txt", package = "adsimR")
  results_adsimR <- compute_selection_coefficients(infile = path, recessive = TRUE, outfile = paste0(tempdir(), "/output.txt"))
  path2 <- system.file("extdata","Selection_scenarios_coefficients.txt", package = "adsimR")
  results_perl <- read.table(path2, header = TRUE)
  expect_equal(results_adsimR, results_perl)
})

test_that("dominant is less than recesive",{
  path <- system.file("extdata","Selection_scenarios.txt", package = "adsimR")
  results_recessive <- compute_selection_coefficients(infile = path, outfile = paste0(tempdir(), "/output.txt"), recessive = TRUE)
  results_dominant <- compute_selection_coefficients(infile = path, outfile = paste0(tempdir(), "/output.txt"), recessive = FALSE)
  expect_equal(results_recessive$s > results_dominant$s, rep(TRUE, 15))
})


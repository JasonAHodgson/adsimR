test_that("test valid input", {
  #create basic dataframe
  df <- data.frame(c = c(1, 2, 3), b = c(0.1, 0.2, 0.3), a = c(0.15, 0.25, 0.35))
  #save as text file
  input_file <- paste0(tempfile("test_input.txt"))
  write.table(df, input_file, sep = "\t", row.names = FALSE, quote = FALSE)
  expect_error(compute_selection_coefficients(infile = input_file,
  outfile = paste0(tempdir(), "/output.txt")), "Input file must contain columns: g, q, qx")
})

test_that("check_uniform_parameters_in_signatures_df stops when there are multiple ksizes", {
  sig_df <- read_signature("G36354.sig.gz")
  expect_error(check_uniform_parameters_in_signatures_df(sig_df))
})

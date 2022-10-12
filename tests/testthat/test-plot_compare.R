test_that("make_compare_mds outputs a data frame of the correct shape", {
  mds <- make_compare_mds(read_compare_csv("tests/testthat/comp_k31.csv"))
  expect_equal(ncol(mds), 3)
  expect_equal(nrow(mds), 6)
})

test_that("make_compare_mds calculates the correct position for one sample", {
  mds <- make_compare_mds(read_compare_csv("tests/testthat/comp_k31.csv"))
  expect_equal(mds[1, 2], -0.37079668)
  expect_equal(mds[1, 3], -0.26357615)
})

test_that("from_list_to_upset_df outputs a correctly formated upset data frame", {
  sample1 <- c("a", "b", "c")
  sample2 <- c("b", "c", "d")
  upset_list <- list(sample1 = sample1, sample2 = sample2)
  upset_df <- from_list_to_upset_df(upset_list)
  expect_equal(as.numeric(rowSums(upset_df)), c(1, 2, 2, 1)) # as.numeric() removes the rownames from rowSums
  expect_equal(colnames(upset_df), c("sample1", "sample2"))
  expect_equal(rownames(upset_df), c("a", "b", "c", "d"))
})

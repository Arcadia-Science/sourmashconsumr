test_that("from_list_to_upset_df outputs a correctly formated upset data frame", {
  sample1 <- c("a", "b", "c")
  sample2 <- c("b", "c", "d")
  upset_list <- list(sample1 = sample1, sample2 = sample2)
  upset_df <- from_list_to_upset_df(upset_list)
  expect_equal(as.numeric(rowSums(upset_df)), c(1, 2, 2, 1)) # as.numeric() removes the rownames from rowSums
  expect_equal(colnames(upset_df), c("sample1", "sample2"))
  expect_equal(rownames(upset_df), c("a", "b", "c", "d"))
})

test_that("from_upset_df_to_intersections records intersections correctly", {
  sample1 <- c("a", "b", "c")
  sample2 <- c("b", "c", "d")
  upset_list <- list(sample1 = sample1, sample2 = sample2)
  upset_df <- from_list_to_upset_df(upset_list)
  intersection_df <- from_upset_df_to_intersections(upset_df)
  # make sure output df has column names
  expect_equal(rownames(intersection_df), c("a", "b", "c", "d"))
  # make sure intersections are accurate
  expect_equal(intersection_df$intersection, c("sample1_0", "sample1_sample2", "sample1_sample2", "0_sample2"))
})

test_that("from_upset_df_to_intersection_summary summarizes correctly", {
  sample1 <- c("a", "b", "c")
  sample2 <- c("b", "c", "d")
  upset_list <- list(sample1 = sample1, sample2 = sample2)
  upset_df <- from_list_to_upset_df(upset_list)
  summary_df <- from_upset_df_to_intersection_summary(upset_df)
  expect_equal(summary_df$n, c(2, 1, 1))
})

test_that("from_upset_df_to_intersection_members pulls out correct intersections", {
  sample1 <- c("a", "b", "c")
  sample2 <- c("b", "c", "d")
  upset_list <- list(sample1 = sample1, sample2 = sample2)
  upset_df <- from_list_to_upset_df(upset_list)
  lst <- from_upset_df_to_intersection_members(upset_df)
  # to pull out values from the list, index into the list using any of the intersection names.
  expect_equal(lst[["sample1_sample2"]], c("b", "c"))
  expect_equal(lst[["sample1_0"]], "a")
})

test_that("make_compare_mds outputs a data frame of the correct shape", {
  mds <- make_compare_mds(read_compare_csv("comp_k31.csv"))
  expect_equal(ncol(mds), 3)
  expect_equal(nrow(mds), 6)
})

test_that("make_compare_mds calculates the correct position for one sample", {
  mds <- make_compare_mds(read_compare_csv("comp_k31.csv"))
  expect_equal(mds[1, 2], -0.37079668)
  expect_equal(mds[1, 3], -0.26357615)
})

test_that("plot_compare_mds produces an object of class ggplot with correct data", {
  mds <- make_compare_mds(read_compare_csv("comp_k31.csv"))
  plt <- plot_compare_mds(mds)
  expect_equal(class(plt)[2], "ggplot")
  expect_equal(plt$data, mds)
})

test_that("plot_compare_mds label boolean behaves correctly by adding second layer", {
  mds <- make_compare_mds(read_compare_csv("comp_k31.csv"))
  plt_true <- plot_compare_mds(mds)
  expect_equal(length(plt_true$layers), 2)
  plt_false <- plot_compare_mds(mds, label = FALSE)
  expect_equal(length(plt_false$layers), 1)
})

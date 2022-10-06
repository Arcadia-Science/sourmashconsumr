
# test read_compare_csv --------------------------------------------------------

test_that("check read_compare_csv is symmetric", {
  df = read_compare_csv("comp_k31.csv", sample_to_rownames = T)
  expect_equal(ncol(df), nrow(df))
})

test_that("check read_csv pass through arguments work correctly for read_compare_csv.", {
  # show_col_types produces a message about info in the tibble.
  expect_message(read_compare_csv("comp_k31.csv", sample_to_rownames = F, show_col_types = T))
})

test_that("check that column names are the same as rownames for read_compare_csv", {
  df = read_compare_csv("comp_k31.csv", sample_to_rownames = T)
  expect_equal(colnames(df), rownames(df))
})


# test read_taxonomy_annotate ---------------------------------------------

test_that("check that read_taxonomy_annotate reads single files genbank db", {
  df_many <- read_taxonomy_annotate(Sys.glob("*genbank*lineages-head*.csv"))
  expect_equal(nrow(df_many), 177)
})

test_that("check that read_taxonomy_annotate reads many files genbank db", {
  df_one <- read_taxonomy_annotate("SRR19888423ass-vs-genbank-2022.03-k31.with-lineages-head23.csv")
  expect_equal(nrow(df_one), 22)
})

test_that("check that read_taxonomy_annotate reads single files gtdb reps db", {
  df_many <- read_taxonomy_annotate(Sys.glob("*gtdbrs207_reps*lineages*.csv"))
  expect_equal(nrow(df_many), 1062)
})

test_that("check that read_taxonomy_annotate reads many files gtdb reps db", {
  df_one <- read_taxonomy_annotate("SRR5947006_gather_gtdbrs207_reps.with-lineages.csv")
  expect_equal(nrow(df_one), 185)
})

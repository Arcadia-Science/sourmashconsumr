
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


# test read_signature -----------------------------------------------------

test_that("read signature works with a signature with only one ksize", {
  # This signature is k=21,abund
  sig_df <- read_signature(file = "reads-s10-s11.sig")
  expect_equal(nrow(sig_df), 770)
})

test_that("read signature works when num is specified", {
  # this signature is num=500,k=21,k=30
  sig_df <- read_signature("genome-s10.fa.gz.sig")
  expect_equal(nrow(sig_df), 2000)
  expect_equal(unique(sig_df$num), 500)
  expect_equal(unique(sig_df$max_hash), 0)
})

test_that("read signature works when there are no abundances", {
  # this signature is no abundances, k21,31,51
  sig_df <- read_signature("TARA_ANW_MAG_00005.sig")
  expect_warning(sig_df$abundances) # should not have an abundance column, so will give warning about uninitialized column
 })

test_that("read signature works when the signature is gzipped", {
  sig_df <- read_signature("G36354.sig.gz")
  expect_equal(nrow(sig_df), 27219)
})

test_that("get_scaled_from_max_hash returns the correct scaled value", {
  sig_df <- read_signature("SRR18071810.sig")
  scaled <- get_scaled_for_max_hash(unique(sig_df$max_hash))
  expect_equal(scaled, 100000)
})

test_that("sigs read with read_signature from diff ver. sourmash can be combined when compliant = TRUE", {
   sig_df_old_sourmash <- read_signature("TARA_ANW_MAG_00005.sig", compliant = T)
   sig_df_new_sourmash <- read_signature("SRR18071810.sig", compliant = T)
   df <- dplyr::bind_rows(sig_df_old_sourmash, sig_df_new_sourmash)
   expect_equal(length(unique(df$filename)), 2)
})


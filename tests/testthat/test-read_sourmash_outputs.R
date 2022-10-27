
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
  # we expect warning messages here given that some of the columns will be missing.
  # because multiple files are read in, multiple warning messages are output
  # expect_warning() only catches one warning
  # suppressWarnings is the recommended way to handle multiple warning you don't care about:
  # https://testthat.r-lib.org/reference/expect_error.html
  suppressWarnings(df_many <- read_taxonomy_annotate(Sys.glob("*genbank*lineages-head*.csv")))
  expect_equal(nrow(df_many), 189)
})

test_that("check that read_taxonomy_annotate reads many files genbank db", {
  expect_warning(
    df_one <- read_taxonomy_annotate("SRR19888423ass-vs-genbank-2022.03-k31.with-lineages-head23.csv"),
    regexp = "n_unique_weighted_found, sum_weighted_found, total_weighted_hashes"
  )
  expect_equal(nrow(df_one), 22)
})

test_that("check that read_taxonomy_annotate reads single files gtdb reps db", {
  # we expect warning messages here given that some of the columns will be missing.
  # because multiple files are read in, multiple warning messages are output
  # expect_warning() only catches one warning
  # suppressWarnings is the recommended way to handle multiple warning you don't care about:
  # https://testthat.r-lib.org/reference/expect_error.html
  suppressWarnings(df_many <- read_taxonomy_annotate(Sys.glob("*gtdbrs207_reps*lineages*.csv")))
  expect_equal(nrow(df_many), 1062)
})

test_that("check that read_taxonomy_annotate reads many files gtdb reps db", {
  expect_warning(
    df_one <- read_taxonomy_annotate("SRR5947006_gather_gtdbrs207_reps.with-lineages.csv"),
    regexp = "n_unique_weighted_found, sum_weighted_found, total_weighted_hashes"
  )
  expect_equal(nrow(df_one), 185)
})


# test read_signature -----------------------------------------------------

test_that("read signature works with a signature with only one ksize", {
  # This signature is k=21,abund
  sig_df <- read_signature(file = "reads-s10-s11.sig")
  expect_equal(nrow(sig_df), 770)
  # make sure we include scaled in the output
  expect_equal(unique(sig_df$scaled), 1000)
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


# sourmash gather ---------------------------------------------------------

test_that("read_gather works with outputs from different versions of sourmash", {
  gather_files <- Sys.glob("*gather.csv")
  # we expect warning messages here given that some of the columns will be missing.
  # because multiple files are read in, multiple warning messages are output
  # expect_warning() only catches one warning
  # suppressWarnings is the recommended way to handle multiple warning you don't care about:
  # https://testthat.r-lib.org/reference/expect_error.html
  suppressWarnings(gather_df <- read_gather(gather_files, intersect_bp_threshold = 0))
  # expect df of correct dimensions when read in with purrr
  expect_equal(ncol(gather_df), 33)
  # lots of columns are missing in this one, so we expect to be warned about it.
  expect_warning(read_gather("47+63_x_gtdb-rs202.gather.csv"),
                 regexp = "The following named parsers don't match the column names: ksize, moltype, scaled, query_n_hashes, query_abundance, query_containment_ani, match_containment_ani, average_containment_ani, max_containment_ani, potential_false_negative, n_unique_weighted_found, sum_weighted_found, total_weighted_hashes"
  )
  # no columns are missing in this one, so we don't expect any problems
  gather_df1 <- read_gather("test1.v450.gather.csv", intersect_bp_threshold = 50000)
  # make sure filtering worked
  expect_true(all(gather_df1$intersect_bp > 50000))
})

test_that("test that specifying no taxonomy_annotate_tibble or file fails with exit code", {
  expect_error(from_taxonomy_annotate_to_metacoder(), regexp = "Neither taxonomy_annotate_df or file were specified. Please specify either taxonomy_annotate_df or file and retry.")
})

test_that("test that user can specify groups with a data frame", {
  run_accessions <- c("SRR5936131", "SRR5947006", "SRR5935765",
                      "SRR5936197", "SRR5946923", "SRR5946920")
  groups <- c("cd", "cd", "cd", "nonibd", "nonibd", "nonibd")
  groups_df <- data.frame(run_accessions, groups)
  expect_message(from_taxonomy_annotate_to_metacoder(file = Sys.glob("*gtdbrs207_reps.with-lineages.csv"),
                                                     groups = groups_df),
                 regexp = "Calculating number of samples with a value greater than 0 for 6 columns in 2 groups for 701 observations")
  expect_message(from_taxonomy_annotate_to_metacoder(file = Sys.glob("*gtdbrs207_reps.with-lineages.csv"),
                                                     groups = groups_df,
                                                     tax_glom_level = "order"),
                 regexp = "Calculating number of samples with a value greater than 0 for 6 columns in 2 groups for 43 observations")

})

test_that("test that to metacoder works with genbank database", {
  expect_message(from_taxonomy_annotate_to_metacoder(file = Sys.glob("*genbank*.with-lineages*.csv")), regexp = "Calculating number of samples with a value greater than 0 for 4 columns for 162 observations")
})

test_that("test that to metacoder fails when an unrecognized level of taxonomy is used", {
  expect_error(from_taxonomy_annotate_to_metacoder(file = Sys.glob("*genbank*.with-lineages*.csv"),
                                                   tax_glom_level = "kingdom"),
                 regexp = "Unrecognized string passed to tax_glom_level. Please use one of species, genus, family, order, class, phylum, or domain.")
})


test_that("test that to metacoder works when given a valid taxonomic level", {
  expect_message(from_taxonomy_annotate_to_metacoder(file = Sys.glob("*genbank*.with-lineages*.csv"),
                                                     tax_glom_level = "order"),
                 regexp = "Summing per-taxon counts from 4 columns for 27 taxa")
})

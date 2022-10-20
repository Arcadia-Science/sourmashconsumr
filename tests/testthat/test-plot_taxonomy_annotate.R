test_that("multiplication works", {
  expect_equal(2 * 2, 4)
})

# taxonomy_annotate_df <- read_taxonomy_annotate(file = Sys.glob("tests/testthat/*gtdbrs207_reps.with-lineages.csv"), separate_lineage = T)
#upset_inputs <- from_taxonomy_annotate_to_upset_inputs(taxonomy_annotate_df, tax_glom_level = "order")
#plot_taxonomy_annotate_upset(upset_inputs, fill = "phylum")

# sapply(rownames(upset_df), function(x){grepl(pattern = x, x = unique(taxonomy_annotate_df$lineage))})
# all(grepl(pattern = rownames(update_df), taxonomy_annotate_df$lineage))
#
# sum(stringr::str_detect(string = unique(taxonomy_annotate_df$lineage), pattern =rownames(upset_df)[1]))
# tmp<- sapply(rownames(upset_df), function(x) {sum(stringr::str_detect(string = unique(taxonomy_annotate_df$lineage), pattern = x))})

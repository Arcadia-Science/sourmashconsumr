
# agglomeration -----------------------------------------------------------

test_that("taxonomic agglomeration works correctly", {
  tax_annot_df1 <- read_taxonomy_annotate(file = Sys.glob("*gtdbrs207_reps.with-lineages.csv"), separate_lineage = T)
  # when no tax glom level is defined, the same df is returned
  tax_glom_df1 <- tax_glom_taxonomy_annotate(tax_annot_df1)
  expect_equal(tax_annot_df1, tax_glom_df1)
  # that it doesn't matter whether tax_annot_df had lineages separated or not
  tax_annot_df2 <- read_taxonomy_annotate(file = Sys.glob("*gtdbrs207_reps.with-lineages.csv"), separate_lineage = F)
  tax_glom_df2 <- tax_glom_taxonomy_annotate(tax_annot_df2)
  expect_equal(tax_annot_df2, tax_glom_df2)
  # tax glom gloms
  tax_glom_df3 <- tax_glom_taxonomy_annotate(tax_annot_df2, tax_glom_level = "phylum")
  expect_equal(length(unique(tax_glom_df3$lineage)), 8)
  # only two levels are in output lineage
  # only run if stringr is installed
  if(rlang::is_installed("stringr") == T){
    semicoloncount <- unique(stringr::str_count(pattern = ";", string = tax_glom_df3$lineage))
    expect_equal(semicoloncount, 1) # at phylum level there should only be 1 semicolon per lineage
  }
})

test_that("make_agglom_cols returns the right columns", {
  expect_equal(make_agglom_cols(tax_glom_level = "phylum", with_query_name = T),
               c("query_name", "domain", "phylum"))
  expect_equal(make_agglom_cols(tax_glom_level = "order", with_query_name = F),
               c("domain", "phylum", "class", "order"))
})
# upset -------------------------------------------------------------------


test_that("from_taxonomy_annotate_to_upset_inputs works correctly", {
  tax_annot_df1 <- read_taxonomy_annotate(file = Sys.glob("*gtdbrs207_reps.with-lineages.csv"), separate_lineage = T)
  upset_inputs <- from_taxonomy_annotate_to_upset_inputs(tax_annot_df1)
  # first two objects should be data frames
  expect_equal(class(upset_inputs[[1]]), "data.frame")
  expect_equal(class(upset_inputs[[2]])[3], "data.frame")
  # tax glom level wasn't set so should be null
  expect_equal(class(upset_inputs[[3]]), "NULL")
  # when no tax glom, all of the rownames in the upset df are lineages in the tax annot df
  expect_true(all(rownames(upset_inputs[[1]]) %in% tax_annot_df1$lineage))
  # when tax glom is specified, all of the rownames in the upset df match the lineages in the output df
  upset_inputs2 <- from_taxonomy_annotate_to_upset_inputs(tax_annot_df1, tax_glom_level = "order")
  expect_true(all(rownames(upset_inputs2[[1]]) %in% upset_inputs2[[2]]$lineage))
  expect_equal(upset_inputs2[[3]], "order")
})

test_that("plot_taxonomy_annotate_upset returns a plot", {
  tax_annot_df1 <- read_taxonomy_annotate(file = Sys.glob("*gtdbrs207_reps.with-lineages.csv"), separate_lineage = T)
  upset_inputs <- from_taxonomy_annotate_to_upset_inputs(tax_annot_df1, "class")
  plt <- plot_taxonomy_annotate_upset(upset_inputs, fill = "phylum")
  expect_equal(class(plt)[[3]], "ggplot")
})


# sankey plot -------------------------------------------------------------

test_that("plot_taxonomy_annotate_sankey returns a plot", {
  taxonomy_annotate_df <- read_taxonomy_annotate(Sys.glob("SRR19*lineage*.csv"), separate_lineage = T)
  # test data has NAs in strain, so this should fail
  expect_error(plot_taxonomy_annotate_sankey(taxonomy_annotate_df = taxonomy_annotate_df),
               regexp = "Some lineages are missing strain information.")
  # but a plot should be returned if a tax glom level is annotated
  sankey_plt <- plot_taxonomy_annotate_sankey(taxonomy_annotate_df = taxonomy_annotate_df,
                                              tax_glom_level = "order")
  expect_equal(class(sankey_plt)[2], "ggplot")
  expect_equal(length(sankey_plt$layers), 3)
  sankey_plt2 <- plot_taxonomy_annotate_sankey(taxonomy_annotate_df = taxonomy_annotate_df,
                                               tax_glom_level = "order",
                                               palette = "#444444")
  expect_equal(length(sankey_plt2$layers), 3)
})

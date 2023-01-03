
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
  # tax glom works for f_unique_to_query too
  tax_glom_df4 <- tax_glom_taxonomy_annotate(tax_annot_df2, tax_glom_level = "phylum", glom_var = "f_unique_to_query")
  # the column names should contain f_unique_to_query
  expect_equal(colnames(tax_glom_df4), c("lineage", "query_name", "f_unique_to_query"))
  # the total sum of f_unique_to_query shouldn't change after tax glom
  expect_equal(sum(tax_annot_df2$f_unique_to_query), sum(tax_glom_df4$f_unique_to_query))
  # there should be fewer rows after agglomeration
  expect_true(nrow(tax_annot_df2) > nrow(tax_glom_df4))
})

test_that("taxnomic agglomeration fails with a variable that isn't an option", {
  tax_annot_df1 <- read_taxonomy_annotate(file = Sys.glob("*gtdbrs207_reps.with-lineages.csv"), separate_lineage = T)
  # when no tax glom level is defined, the same df is returned
  expect_error(tax_glom_taxonomy_annotate(tax_annot_df1, tax_glom_level = "order", glom_var = "std_abund"),
               regexp = "The variable you supplied is not a valid glom_var.")
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
  # check that it still returns a plot when no tax glom was performed initially
  upset_inputs2 <- from_taxonomy_annotate_to_upset_inputs(tax_annot_df1)
  plt2 <- plot_taxonomy_annotate_upset(upset_inputs2, fill = NULL)
  expect_equal(class(plt2)[[3]], "ggplot")

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

# ts alluvial -------------------------------------------------------------

test_that("plot_taxonomy_annotate_ts_alluvial returns a plot", {
  taxonomy_annotate_df <- read_taxonomy_annotate(Sys.glob("SRR19*lineage*.csv"), separate_lineage = T)
  time_df <- data.frame(query_name = unique(taxonomy_annotate_df$query_name),
                        time = c(1, 2, 3, 4))
  alluvial_plt <- plot_taxonomy_annotate_ts_alluvial(taxonomy_annotate_df, time_df = time_df, tax_glom_level = "genus")
  expect_equal(class(alluvial_plt)[2], "ggplot")
  expect_equal(length(alluvial_plt$layers), 2)
})

# multi strain detection --------------------------------------------------------

test_that("from_taxonomy_annotate_to_multi_strains() returns a plot", {
  taxonomy_annotate_df <- read_taxonomy_annotate(Sys.glob("SRR19*lineage*.csv"), separate_lineage = T)
  lst <- from_taxonomy_annotate_to_multi_strains(taxonomy_annotate_df = taxonomy_annotate_df)
  # the mapping should be to f_match if plotting worked
  expect_equal(c("$", ".data", "f_match"), as.character(rlang::quo_get_expr(lst$plt$layers[[1]]$mapping$size)))
  # only one query (SRR19888434) and one species (Brevibacterium aurantiacum) have an f_match over 1.1
  expect_equal(1, nrow(lst$candidate_species_with_multiple_strains))
  # there were six genomes detected
  expect_equal(6, nrow(lst$plt_data))
  # test that the stop message is returned when there are no f_match results
  taxonomy_annotate_df <- read_taxonomy_annotate(Sys.glob("SRR19888427*lineage*.csv"), separate_lineage = T)
  lst <- from_taxonomy_annotate_to_multi_strains(taxonomy_annotate_df = taxonomy_annotate_df)
  expect_equal(0, nrow(lst$candidate_species_with_multiple_strains))
})

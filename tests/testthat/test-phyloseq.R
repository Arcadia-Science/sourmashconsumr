test_that("from_taxonomy_annotate_to_tax_table produces a taxonomy table", {
  tax_annot_gtdb <- read_taxonomy_annotate(file = Sys.glob("*gtdbrs207_reps.with-lineages.csv"), separate_lineage = T)
  tax_table_gtdb <- from_taxonomy_annotate_to_tax_table(tax_annot_gtdb)
  expect_equal(colnames(tax_table_gtdb), c('domain', 'phylum', 'class', 'order', 'family', 'genus', 'species', 'strain'))
  expect_equal(nrow(tax_table_gtdb), 437)
  expect_true(all(is.na(tax_table_gtdb$strain))) # GTDB doesn't have strain information
  tax_annot_genbank <- read_taxonomy_annotate(Sys.glob("*genbank*lineages-head*.csv"))
  tax_table_genbank <- from_taxonomy_annotate_to_tax_table(tax_annot_genbank)
  expect_equal(colnames(tax_table_genbank), c('domain', 'phylum', 'class', 'order', 'family', 'genus', 'species', 'strain'))
  expect_equal(nrow(tax_table_genbank), 110)
})

test_that("from_taxonomy_annotate_to_count_table produces a count table", {
  tax_annot_gtdb <- read_taxonomy_annotate(file = Sys.glob("*gtdbrs207_reps.with-lineages.csv"), separate_lineage = T)
  count_table_gtdb <- from_taxonomy_annotate_to_count_table(tax_annot_gtdb)
  expect_equal(colSums(count_table_gtdb),
               c('SRR5935765' = 952790, 'SRR5936131' = 2352470,
                 'SRR5936197' = 2083110, 'SRR5946920' = 535214,
                 'SRR5946923' = 1281110, 'SRR5947006' = 1203502))
  tax_annot_genbank <- read_taxonomy_annotate(Sys.glob("*genbank*lineages-head*.csv"))
  count_table_genbank <- from_taxonomy_annotate_to_count_table(tax_annot_genbank)
  expect_equal(colSums(count_table_genbank),
               c('SRR19888423' = 25012, 'SRR19888427' = 32318,
                 'SRR19888432' = 27645, 'SRR19888434' = 31283))
})

test_that("from_taxonomy_annotate_to_phyloseq", {
  # check for gtdb
  tax_annot_gtdb <- read_taxonomy_annotate(file = Sys.glob("*gtdbrs207_reps.with-lineages.csv"), separate_lineage = T)
  metadata <- data.frame(group = c(rep("group1", 3), rep("group2", 3)))
  rownames(metadata) <- unique(tax_annot_gtdb$query_name)
  physeq_gtdb <- from_taxonomy_annotate_to_phyloseq(taxonomy_annotate_df = tax_annot_gtdb, metadata_df = metadata)
  expect_equal(class(physeq_gtdb)[1], "phyloseq")
  # check for genbank
  tax_annot_genbank <- read_taxonomy_annotate(file = Sys.glob("*genbank*lineages-head*.csv"), separate_lineage = T)
  metadata <- data.frame(group = c(rep("group1", 2), rep("group2", 2)))
  rownames(metadata) <- unique(tax_annot_genbank$query_name)
  physeq_genbank <- from_taxonomy_annotate_to_phyloseq(taxonomy_annotate_df = tax_annot_genbank, metadata_df = metadata)
  expect_equal(class(physeq_genbank)[1], "phyloseq")
  # try with no metadata
  physeq_genbank2 <- from_taxonomy_annotate_to_phyloseq(taxonomy_annotate_df = tax_annot_genbank)
  expect_equal(class(physeq_genbank2)[1], "phyloseq")
})

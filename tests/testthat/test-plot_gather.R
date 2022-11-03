test_that("plot_gather_classified produces a barchart", {
  gather_files <- Sys.glob("*gather.csv")
  gather_df <- read_gather(gather_files, intersect_bp_threshold = 0)
  bar_plt <- plot_gather_classified(gather_df)
  expect_equal(class(bar_plt), c("gg", "ggplot"))

})

test_that("from_gather_to_upset_df creates a properly formatted upset df", {
  gather_files <- Sys.glob("*gather.csv")
  gather_df <- read_gather(gather_files, intersect_bp_threshold = 0)
  upset_df <- from_gather_to_upset_df(gather_df)
  expect_equal(nrow(upset_df), length(unique(gather_df$genome_accession))) # each genome should become a row
  expect_equal(colnames(upset_df), c("47+63", "MAG3_1", "SRR19888434", "SRR19888438", "SRR19888440", "test1"))
})

test_that("plot_gather_upset produces a complexupset plot", {
  gather_files <- Sys.glob("*gather.csv")
  gather_df <- read_gather(gather_files, intersect_bp_threshold = 0)
  upset_df <- from_gather_to_upset_df(gather_df)
  expect_equal(class(plot_gather_upset(upset_df)), c("patchwork", "gg", "ggplot"))
  # this should fail with this gather_df because we repeat some genomes across checkV databases
  expect_error(plot_gather_upset(upset_df, color_by_database = T, gather_df = gather_df), regexp = "At least one genome accession was in multiple databases for the supplied gather results. Either use color_by_database == FALSE or limit which samples you plot.")
  # if we limit which gather files we read in, we should be able to get a plot
  gather_files <- c("SRR19888440ass-vs-genbank-2022.03-k31-gather.csv",
                    "lemonade-MAG3.x.gtdb.gather.csv",
                    "test1.v450.gather.csv")
  gather_df <- read_gather(gather_files, intersect_bp_threshold = 0)
  upset_df <- from_gather_to_upset_df(gather_df)
  upset_plt <- plot_gather_upset(upset_df, color_by_database = TRUE, gather_df = gather_df )
  expect_equal(class(upset_plt), c("patchwork", "gg", "ggplot"))
})

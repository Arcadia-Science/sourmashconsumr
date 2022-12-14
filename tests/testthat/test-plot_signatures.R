
# signatures to upset plot ------------------------------------------------

test_that("check_uniform_parameters_in_signatures_df stops when there are multiple ksizes", {
  sig_df <- read_signature("G36354.sig.gz")
  expect_error(check_uniform_parameters_in_signatures_df(sig_df))
})

test_that("from_signatures_to_upset_df produces an upset df", {
  sig_df <- Sys.glob("SRR*ass.sig*") %>%
    purrr::map_dfr(read_signature)
  # should fail because we have multiple k-mer sizes
  expect_error(from_signatures_to_upset_df(sig_df))
  sig_df <- sig_df %>%
    dplyr::filter(ksize == 21)
  upset_df <- from_signatures_to_upset_df(sig_df)
  expect_equal(nrow(upset_df), length(unique(sig_df$mins))) # each minhash should become a row
  expect_equal(colnames(upset_df), c('SRR19888423', 'SRR19888427', 'SRR19888432', 'SRR19888434', 'SRR19888438', 'SRR19888440'))
})

test_that("from_signatures_to_upset_df produces an upset df", {
  sig_df <- Sys.glob("SRR*ass.sig*") %>%
    purrr::map_dfr(read_signature) %>%
    dplyr::filter(ksize == 21)
  upset_df <- from_signatures_to_upset_df(sig_df)
  plt <- plot_signatures_upset(upset_df)
  expect_equal(length(plt$layers), 4)
})

# signatures to rarefaction curves ----------------------------------------

test_that("from_signatures_to_rarefaction_df runs", {
  sig_df <- read_signature(Sys.glob("SRR*ass.sig*"))
  # should fail because we have multiple k-mer sizes
  expect_error(from_signatures_to_rarefaction_df(sig_df))
  sig_df <- sig_df %>%
    dplyr::filter(ksize == 21)
  rarecurve_df <- from_signatures_to_rarefaction_df(sig_df, step = 1)
  expect_equal(max(rarecurve_df$num_kmers_sampled), 684)
  expect_equal(max(rarecurve_df$num_kmers_observed), 637)
})

test_that("plot_signatures_rarefaction produces ggplot with right info", {
  sig_df <- read_signature(Sys.glob("SRR*ass.sig*")) %>%
    dplyr::filter(ksize == 21)
  rarecurve_df <- from_signatures_to_rarefaction_df(sig_df, step = 1)
  plt1 <- plot_signatures_rarefaction(rarecurve_df, fraction_of_points_to_plot = 1)
  plt2 <- plot_signatures_rarefaction(rarecurve_df, fraction_of_points_to_plot = 2)
  expect(nrow(plt1$data) > nrow(plt2$data), ok = T) # the data underlying the first plot should have more rows than the second
  expect_equal(plt1$theme, plt2$theme) # theme should match between plots
  plt1 <- plt1 + ggplot2::theme_minimal()
  expect(plt1$theme$panel.background != plt2$theme$panel.background, ok = T) # after changing a ggplot layer themes should no longer be equal
})

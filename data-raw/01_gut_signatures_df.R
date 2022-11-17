gut_signatures_df <- Sys.glob("inst/extdata/SRR*sig.gz") %>%
  purrr::map_dfr(sourmashconsumr::read_signature)
usethis::use_data(gut_signatures_df, overwrite = TRUE, compress = "xz")

gut_gather_df <- sourmashconsumr::read_gather(Sys.glob("inst/extdata/SRR*reps.csv.gz"))
usethis::use_data(gut_gather_df, overwrite = TRUE, compress = "xz")

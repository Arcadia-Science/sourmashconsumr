gut_compare_df <- sourmashconsumr::read_compare_csv("inst/extdata/gut_compare.csv.gz")
usethis::use_data(gut_compare_df, overwrite = TRUE, compress = "xz")

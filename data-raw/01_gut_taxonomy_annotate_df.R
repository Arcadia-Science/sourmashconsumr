gut_taxonomy_annotate_df <- sourmashconsumr::read_taxonomy_annotate(Sys.glob("inst/extdata/SRR*with-lineages*csv.gz"))
usethis::use_data(gut_taxonomy_annotate_df, overwrite = TRUE, compress = "xz")

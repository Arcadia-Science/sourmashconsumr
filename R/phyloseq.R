#' Transform the output of sourmash taxonomy annotate into a phyloseq tax table
#'
#' @param taxonomy_annotate_df Data frame containing outputs from sourmash taxonomy annotate.
#'
#' @return A phyloseq tax_table object.
#' @export
#'
#' @examples
#' \dontrun{
#' from_taxonomy_annotate_to_tax_table(taxonomy_annotate_df)
#' }
from_taxonomy_annotate_to_tax_table <- function(taxonomy_annotate_df) {
  tax_table <- taxonomy_annotate_df %>%
    dplyr::select("name", "lineage") %>%
    dplyr::distinct() %>%
    tidyr::separate("lineage", into = c("domain", "phylum", "class", "order", "family", "genus", "species", "strain"),
                    sep = ";", extra = "drop", fill = "right") %>%
    tibble::column_to_rownames("name")
  return(tax_table)
}

#' Transform the output of sourmash taxonomy annotate into a phyloseq OTU table (count table)
#'
#' @param taxonomy_annotate_df Data frame containing outputs from sourmash taxonomy annotate.
#'
#' @return A phyloseq otu_table object
#' @export
#'
#' @examples
#' \dontrun{
#' from_taxonomy_annotate_to_count_table(taxonomy_annotate_df)
#' }
from_taxonomy_annotate_to_count_table <- function(taxonomy_annotate_df) {
  count_table <- taxonomy_annotate_df %>%
    dplyr::select("query_name", "name", "n_unique_kmers") %>% # select only the columns that have information we need
    tidyr::pivot_wider(id_cols = "name", names_from = "query_name", values_from = "n_unique_kmers") # transform to wide format

  count_table[is.na(count_table)] <- 0 # replace all NAs with 0

  count_table <- count_table %>%
    tibble::column_to_rownames("name") # move the metagenome sample name to a rowname

  return(count_table)
}

#' Transform the output of sourmash taxonomy annotate into a phyloseq object
#'
#' @description
#' `from_taxonomy_annotate_to_phyloseq()` transforms a data frame that contains sourmash taxonomy annotate results into a phyloseq object.
#' Counts are derived from the `n_unique_kmers` column, the abundance-weighted number of unique k-mers overlapping between a query and its match in the database.
#'
#' @param taxonomy_annotate_df Data frame containing outputs from sourmash taxonomy annotate.
#' @param metadata_df Optional. A data frame of metadata. The row names must match the values of `query_names` in the `taxonomy_annotate_df`. Any additional columns may be added.
#'
#' @return A phyloseq object.
#' @export
#'
#' @examples
#' \dontrun{
#' from_taxonomy_annotate_to_phyloseq(taxonomy_annotate_df)
#' }
from_taxonomy_annotate_to_phyloseq <- function(taxonomy_annotate_df, metadata_df = NULL) {
  # make tax_table and count_table objects
  tax_table <- from_taxonomy_annotate_to_tax_table(taxonomy_annotate_df)
  count_table <- from_taxonomy_annotate_to_count_table(taxonomy_annotate_df)

  # build phyloseq object
  if(!is.null(metadata_df)){
    # add metadata checks
    if(!all(rownames(metadata_df) %in% taxonomy_annotate_df$query_name)){
      stop("Not all metadata_df samples (row names) occur in the taxonomy_annotate_df (query_name).")
    }
    if(!all(taxonomy_annotate_df$query_name %in% rownames(metadata_df))){
      stop("Not all samples in taxonomy_annotate_df (query_name) have metadata_df samples (rown ames).")
    }

    phyloseq_obj <- phyloseq::phyloseq(phyloseq::otu_table(count_table, taxa_are_rows = T),
                                       phyloseq::tax_table(as.matrix(tax_table)),
                                       phyloseq::sample_data(metadata_df))
  } else if(is.null(metadata_df)){
    phyloseq_obj <- phyloseq::phyloseq(phyloseq::otu_table(count_table, taxa_are_rows = T),
                                       phyloseq::tax_table(as.matrix(tax_table)))
  }
  return(phyloseq_obj)
}

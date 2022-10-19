#' Agglomerate counts of same lineage to specified level of taxonomy.
#'
#' @description
#' Inspired by phyloseq::tax_glom(), this method summarizes k-mer counts from genomes that have the same taxonomy at a user-specified taxonomy rank.
#' Agglomeration occurs within each sample, meaning the number of k-mers is only summed within each query_name.
#' This function returns a data frame with the columns lineage, query_name, and n_unique_kmers.
#'
#' @param taxonomy_annotate_df Data frame containing outputs from sourmash taxonomy annotate. Can contain results from one or many runs of sourmash taxonomy annotate. Agglomeration occurs within each query.
#' @param tax_glom_level Character. NULL by default, meaning no agglomeration is done. Valid options are "domain", "phylum", "class", "order", "family", "genus", and "species". When a valid option is supplied, k-mer counts are agglomerated to that level
#'
#' @return A data frame.
#' @export
#'
#' @examples
tax_glom_taxonomy_annotate <- function(taxonomy_annotate_df, tax_glom_level = NULL){
  # if tax_glom_level is not defined, return the taxonomy_annotate_df unchanged
  if(is.null(tax_glom_level)){
    message("tax_glom_level not defined, not doing any agglomeration. Returning input data frame unchanged.")
    return(taxonomy_annotate_df)
  }

  # if tax_glom_level is defined, parse lineage and agglomerate counts to that level of taxonomy
  if(!is.null(tax_glom_level)){
    # make sure only except arguments are returned
    if(!tax_glom_level %in% c("domain", "phylum", "class", "order", "family", "genus", "species")){
      stop("Unrecognized string passed to tax_glom_level. Please use one of species, genus, family, order, class, phylum, or domain.")
    }
    # agglomerate to level of taxonomy
    if(tax_glom_level == "domain"){
      agglom_cols <- c("query_name", "domain")
    } else if(tax_glom_level == "phylum"){
      agglom_cols <- c("query_name", "domain", "phylum")
    } else if(tax_glom_level == "class"){
      agglom_cols <- c("query_name", "domain", "phylum", "class")
    } else if(tax_glom_level == "order"){
      agglom_cols <- c("query_name", "domain", "phylum", "class", "order")
    } else if(tax_glom_level == "family"){
      agglom_cols <- c("query_name", "domain", "phylum", "class", "order", "family")
    } else if(tax_glom_level == "genus"){
      agglom_cols <- c("query_name", "domain", "phylum", "class", "order", "family", "genus")
    } else if(tax_glom_level == "species"){
      agglom_cols <- c("query_name", "domain", "phylum", "class", "order", "family", "genus", "species")
    }
  }

  taxonomy_annotate_df <- taxonomy_annotate_df %>%
    dplyr::select(genome_accession, lineage, query_name, n_unique_kmers) %>%
    tidyr::separate(lineage, into = c("domain", "phylum", "class", "order", "family", "genus", "species", "strain"), sep = ";", remove = F, fill = "right") %>% # split the lineage to each taxonomy rank
    dplyr::group_by_at(dplyr::vars(dplyr::all_of(agglom_cols))) %>% # group using the lineage columns up to the user-specified taxonomic rank
    dplyr::summarize(n_unique_kmers = sum(n_unique_kmers)) %>%
    dplyr::ungroup() %>%
    tidyr::unite(lineage, all_of(agglom_cols[-1]), sep = ";", remove = TRUE) %>% # make a new lineage string that includes only the agglom columns
    dplyr::select(lineage, query_name, n_unique_kmers)

  return(taxonomy_annotate_df)
}

#' Title
#'
#' @description
#' The file parameter was removed because the taxonomy_annotate_df is needed to add fill color to the upset plot.
#' Using the file parameter would bypass the creation of this object in the global environment.
#'
#' @param taxonomy_annotate_df
#'
#'
#' @return
#' @export
#'
#' @examples
from_taxonomy_annotate_to_upset_df <- function(taxonomy_annotate_df,
                                               tax_glom_level = NULL){

  if(!is.null(tax_glom_level)){
    taxonomy_annotate_df <- tax_glom_taxonomy_annotate(taxonomy_annotate_df, tax_glom_level = tax_glom_level)
  }

  # select only columns that are relevant
  taxonomy_annotate_df <- taxonomy_annotate_df %>%
    dplyr::select(lineage, query_name)

  # turn data frame into list.
  # Each index in the list is named after the query_name it represents.
  # Each index contains a vector of lineages that were identified in the sample
  upset_list <- taxonomy_annotate_df %>%
    dplyr::group_by(query_name) %>%
    {setNames(dplyr::group_split(.), dplyr::group_keys(.)[[1]])} # absolutely whack line of code that allows us to set the list names by the group_by variable

  # extract just the minhashes
  upset_list <- lapply(upset_list, function(x) {x$lineage})

  # convert list into upset-compliant data.frame
  upset_df <- from_list_to_upset_df(upset_list)

  return(upset_df)
}

#' Title
#'
#' @description
#' If fill is specified, `taxonomy_annotate_df` must also specified as it is used to derive metadata passed to fill.
#' @param upset_df
#' @param taxonomy_annotate_df
#' @param fill
#'
#' @return
#' @export
#'
#' @examples
plot_taxonomy_annotate_upset <- function(upset_df, taxonomy_annotate_df = NULL, fill = NULL){
  if(!is.null(fill) & is.null(taxonomy_annotate_df)){
    stop("taxonomy_annotate_df is used to derive fill information. Please supply the taxonomy_annotate_df used to make the upset_df.")
  }

  if(!is.null(taxonomy_annotate_df) & is.null(fill)){
    message("taxonomy_annotate_df supplied by fill not specified. That's fine! But you won't see any pretty colors. Add a fill variable to get colors")
  }

  if(!is.null(fill) & !is.null(taxonomy_annotate_df)){
    # make sure that taxonomy_annotate_df has all of the variables in upset_df.
    # If not, the wrong taxonomy_annotate_df was probably supplied
    check_vars_match_upset_df_and_taxonomy_annotate_df <- sapply(rownames(upset_df), function(x) {sum(stringr::str_detect(string = unique(taxonomy_annotate_df$lineage), pattern = x))})
    # stop if every rowname isn't in the original taxonomy_annotate_df
    stopifnot(all(check_vars_match_upset_df_and_taxonomy_annotate_df > 0))
  }

  # plot the upset plot
  upset(upset_df, intersect = names(sourmash_taxonomy_upset_list), set_sizes = F,
        base_annotations=list(
          '# lineages'=intersection_size(text=list(vjust=0.4, hjust=.05, angle=90),
                                         text_colors=c(on_background='black', on_bar='black'),
                                         mapping=aes(fill=database)) +
            scale_fill_brewer(palette = "Set2"))
  )

}


# upset_df <- from_taxonomy_annotate_to_upset_df(file = Sys.glob("tests/testthat/*gtdbrs207_reps.with-lineages.csv"))
# taxonomy_annotate_df <- read_taxonomy_annotate(file = Sys.glob("tests/testthat/*gtdbrs207_reps.with-lineages.csv"))
# sapply(rownames(upset_df), function(x){grepl(pattern = x, x = unique(taxonomy_annotate_df$lineage))})
# all(grepl(pattern = rownames(update_df), taxonomy_annotate_df$lineage))
#
# sum(stringr::str_detect(string = unique(taxonomy_annotate_df$lineage), pattern =rownames(upset_df)[1]))
# tmp<- sapply(rownames(upset_df), function(x) {sum(stringr::str_detect(string = unique(taxonomy_annotate_df$lineage), pattern = x))})


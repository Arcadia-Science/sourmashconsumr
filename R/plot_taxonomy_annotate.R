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
from_taxonomy_annotate_to_upset_inputs <- function(taxonomy_annotate_df,
                                                   tax_glom_level = NULL){

  if(!is.null(tax_glom_level)){
    taxonomy_annotate_df <- tax_glom_taxonomy_annotate(taxonomy_annotate_df, tax_glom_level = tax_glom_level)
  }

  # select only columns that are relevant
  taxonomy_annotate_df <- taxonomy_annotate_df %>%
    dplyr::select(lineage, query_name) %>%
    dplyr::distinct()

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

  # combine the processed taxonomy_annotate_df with the upset_df so that the processed taxonomy_annotate_df can be used to add color to the upset plot
  upset_inputs <- list(upset_df = upset_df, taxonomy_annotate_df = taxonomy_annotate_df, tax_glom_level = tax_glom_level)

  return(upset_inputs)
}

#' Title
#'
#' @description

#' @param upset_inputs
#' @param fill
#'
#' @return
#' @export
#'
#' @examples
plot_taxonomy_annotate_upset <- function(upset_inputs, fill = NULL){
  upset_df <- upset_inputs[[1]]
  taxonomy_annotate_df <- upset_inputs[[2]]
  tax_glom_level <- upset_inputs[[3]]

  if(!is.null(fill)){
    # figure out which columns can encode fill color from taxglom
    if(tax_glom_level == "domain"){
      fill_cols <- c("domain")
    } else if(tax_glom_level == "phylum"){
      fill_cols <- c("domain", "phylum")
    } else if(tax_glom_level == "class"){
      fill_cols <- c("domain", "phylum", "class")
    } else if(tax_glom_level == "order"){
      fill_cols <- c("domain", "phylum", "class", "order")
    } else if(tax_glom_level == "family"){
      fill_cols <- c("domain", "phylum", "class", "order", "family")
    } else if(tax_glom_level == "genus"){
      fill_cols <- c("domain", "phylum", "class", "order", "family", "genus")
    } else if(tax_glom_level == "species"){
      fill_cols <- c("domain", "phylum", "class", "order", "family", "genus", "species")
    }

    # join together upset df with metadata in taxonomy_annotate_df
    taxonomy_annotate_df <-  taxonomy_annotate_df %>%
      dplyr::select(lineage) %>%
      tidyr::separate(lineage, into = fill_cols, sep = ";", remove = F) %>%
      dplyr::distinct() %>%
      dplyr::select(lineage, fill = tidyselect::all_of(fill)) # make new column name based on the identity of parameter fill

    upset_df <- upset_df %>%
      tibble::rownames_to_column("lineage") %>%
      dplyr::left_join(taxonomy_annotate_df, by = "lineage") %>%
      tibble::column_to_rownames("lineage")
  }

  # plot the upset plot
  plt <- ComplexUpset::upset(upset_df, intersect = unique(upset_inputs[[2]]$query_name), set_sizes = F,
                             base_annotations=list(
                               '# lineages'=ComplexUpset::intersection_size(text=list(vjust=0.4, hjust=.05, angle=90),
                                                                            text_colors=c(on_background='black', on_bar='black'),
                                                                            mapping=ggplot2::aes(fill=fill)) +
                                 ggplot2::scale_fill_brewer(palette = "Set2"))
  )
  return(plt)
}

# taxonomy_annotate_df <- read_taxonomy_annotate(file = Sys.glob("tests/testthat/*gtdbrs207_reps.with-lineages.csv"), separate_lineage = T)
#upset_inputs <- from_taxonomy_annotate_to_upset_inputs(taxonomy_annotate_df, tax_glom_level = "order")
# plot_taxonomy_annotate_upset(upset_inputs, fill = "phylum")

# sapply(rownames(upset_df), function(x){grepl(pattern = x, x = unique(taxonomy_annotate_df$lineage))})
# all(grepl(pattern = rownames(update_df), taxonomy_annotate_df$lineage))
#
# sum(stringr::str_detect(string = unique(taxonomy_annotate_df$lineage), pattern =rownames(upset_df)[1]))
# tmp<- sapply(rownames(upset_df), function(x) {sum(stringr::str_detect(string = unique(taxonomy_annotate_df$lineage), pattern = x))})

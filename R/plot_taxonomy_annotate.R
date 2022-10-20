#' Agglomerate counts of same lineage to specified level of taxonomy.
#'
#' @description
#' Inspired by phyloseq::tax_glom(), this method summarizes k-mer counts from genomes that have the same taxonomy at a user-specified taxonomy rank.
#' Agglomeration occurs within each sample, meaning the number of k-mers is only summed within each query_name.
#' This function returns a data frame with the columns `lineage`, `query_name`, and `n_unique_kmers`.
#'
#' @param taxonomy_annotate_df Data frame containing outputs from sourmash taxonomy annotate. Can contain results from one or many runs of sourmash taxonomy annotate. Agglomeration occurs within each query.
#' @param tax_glom_level Character. NULL by default, meaning no agglomeration is done. Valid options are "domain", "phylum", "class", "order", "family", "genus", and "species". When a valid option is supplied, k-mer counts are agglomerated to that level
#'
#' @return A data frame.
#' @export
#'
#' @importFrom rlang .data
#'
#' @examples
#' \dontrun{
#' tax_glom_taxonomy_annotate()
#' }
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
    dplyr::select("genome_accession", "lineage", "query_name", "n_unique_kmers") %>%
    tidyr::separate(.data$lineage, into = c("domain", "phylum", "class", "order", "family", "genus", "species", "strain"), sep = ";", remove = F, fill = "right") %>%
    dplyr::group_by_at(dplyr::vars(dplyr::all_of(agglom_cols))) %>%
    dplyr::summarize(n_unique_kmers = sum(.data$n_unique_kmers)) %>%
    dplyr::ungroup() %>%
    tidyr::unite(col = "lineage", tidyselect::all_of(agglom_cols[-1]), sep = ";", remove = TRUE) %>%
    dplyr::select("lineage", "query_name", "n_unique_kmers")

  return(taxonomy_annotate_df)
}

#' Transform a taxonomy annotate data frame into an upset plot compliant data frame
#'
#' @description
#' `from_taxonomy_annotate_to_upset_inputs()` transforms a data frame with produced using read_taxonomy_annotate on many results produced by sourmash taxonomy annotate into a upset-compliant data frame.
#' The function can optionally agglomerate to different levels of taxonomic rank (e.g. phylum) and then produce the upset compliant data frame after agglomeration.
#'
#' @param taxonomy_annotate_df Data frame containing outputs from sourmash taxonomy annotate. Can contain results from one or many runs of sourmash taxonomy annotate. Agglomeration occurs within each query.
#' @param tax_glom_level Optional character string specifying the taxonomic rank to agglomerate k-mer counts. Must be one of "domain", "phylum", "class", "order", "family", "genus", "species."
#'
#' @return A list. The first object in the list is an upset compliant data frame. The second is the taxonomy_annotate_df used to build the upset data frame. The last is the tax_glom_level.
#' @export
#'
#' @importFrom rlang .data
#'
#' @examples
#' \dontrun{
#' from_taxonomy_annotate_to_upset_inputs()
#' }
from_taxonomy_annotate_to_upset_inputs <- function(taxonomy_annotate_df,
                                                   tax_glom_level = NULL){

  if(!is.null(tax_glom_level)){
    taxonomy_annotate_df <- tax_glom_taxonomy_annotate(taxonomy_annotate_df, tax_glom_level = tax_glom_level)
  }

  # select only columns that are relevant
  taxonomy_annotate_df <- taxonomy_annotate_df %>%
    dplyr::select("lineage", "query_name") %>%
    dplyr::distinct()

  # turn data frame into list.
  # Each index in the list is named after the query_name it represents.
  # Each index contains a vector of lineages that were identified in the sample
  . <- NULL # need to set . as global var because of whack line.
  # I don't actually like this approach, but I think it's ok for now until I can find a workaround for said line.
  upset_list <- taxonomy_annotate_df %>%
    dplyr::group_by(.data$query_name) %>%
    {stats::setNames(dplyr::group_split(.), dplyr::group_keys(.)[[1]])} # absolutely whack line of code that allows us to set the list names by the group_by variable

  # extract just the minhashes
  upset_list <- lapply(upset_list, function(x) {x$lineage})

  # convert list into upset-compliant data.frame
  upset_df <- from_list_to_upset_df(upset_list)

  # combine the processed taxonomy_annotate_df with the upset_df so that the processed taxonomy_annotate_df can be used to add color to the upset plot
  upset_inputs <- list(upset_df = upset_df, taxonomy_annotate_df = taxonomy_annotate_df, tax_glom_level = tax_glom_level)

  return(upset_inputs)
}

#' Plot an upset plot of lineage intersections between samples
#'
#' @description
#' `plot_taxonomy_annotate_upset()` produces a ComplexUpset plot displaying the intersection of taxonomic lineages observed in many samples.
#' The bar chart that displays the intersection size between samples can optionally be colored by taxonomy lineage (e.g. phylum).
#'
#' @param upset_inputs List of inputs produced by from_taxonomy_annotate_to_upset_inputs().
#' @param fill Optional argument specifying which level of taxonomy to fill the upset plot intersections with. Uses the Set2 palette so cannot visualize more than 8 levels.
#'
#' @return A ComplexUpset plot
#' @export
#'
#' @importFrom rlang .data
#'
#' @examples
#' \dontrun{
#' plot_taxonomy_annotate_upset()
#' }
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
      dplyr::select("lineage") %>%
      tidyr::separate("lineage", into = fill_cols, sep = ";", remove = F) %>%
      dplyr::distinct() %>%
      dplyr::select("lineage", fill = tidyselect::all_of(fill)) # make new column name based on the identity of parameter fill

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
                                                                            mapping=ggplot2::aes(fill = .data$fill)) +
                                 ggplot2::scale_fill_brewer(palette = "Set2") +
                                 ggplot2::labs(fill = fill))
  )
  return(plt)
}

# taxonomy agglomeration --------------------------------------------------
#' Helper function to define column names for taxonomy agglomeration
#'
#' @param tax_glom_level Character. Valid options are "domain", "phylum", "class", "order", "family", "genus", and "species".
#' @param with_query_name Boolean indicating whether the column name "query_name" should be included.
#'
#' @return A character vector of column names to be used for agglomeration
#' @export
#'
#' @examples
#' make_agglom_cols("species", with_query_name = FALSE)
make_agglom_cols <- function(tax_glom_level, with_query_name = FALSE){
  if(tax_glom_level == "domain"){
    agglom_cols <- c("domain")
  } else if(tax_glom_level == "phylum"){
    agglom_cols <- c("domain", "phylum")
  } else if(tax_glom_level == "class"){
    agglom_cols <- c("domain", "phylum", "class")
  } else if(tax_glom_level == "order"){
    agglom_cols <- c("domain", "phylum", "class", "order")
  } else if(tax_glom_level == "family"){
    agglom_cols <- c("domain", "phylum", "class", "order", "family")
  } else if(tax_glom_level == "genus"){
    agglom_cols <- c("domain", "phylum", "class", "order", "family", "genus")
  } else if(tax_glom_level == "species"){
    agglom_cols <- c("domain", "phylum", "class", "order", "family", "genus", "species")
  }

  if(with_query_name == TRUE){
    agglom_cols <- c("query_name", agglom_cols)
    return(agglom_cols)
  }
  if(with_query_name == FALSE){
    return(agglom_cols)
  }
}

#' Agglomerate counts of same lineage to specified level of taxonomy.
#'
#' @description
#' Inspired by phyloseq::tax_glom(), this method summarizes k-mer counts from genomes that have the same taxonomy at a user-specified taxonomy rank.
#' Agglomeration occurs within each sample, meaning the number of k-mers is only summed within each query_name.
#' This function returns a data frame with the columns `lineage`, `query_name`, and `n_unique_kmers`.
#'
#' @param taxonomy_annotate_df Data frame containing outputs from sourmash taxonomy annotate. Can contain results from one or many runs of sourmash taxonomy annotate. Agglomeration occurs within each query.
#' @param tax_glom_level Character. NULL by default, meaning no agglomeration is done. Valid options are "domain", "phylum", "class", "order", "family", "genus", and "species". When a valid option is supplied, k-mer counts are agglomerated to that level
#' @param glom_var Character. One of "n_unique_kmers" or "f_unique_to_query".
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
tax_glom_taxonomy_annotate <- function(taxonomy_annotate_df, tax_glom_level = NULL, glom_var = "n_unique_kmers"){
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
    agglom_cols <- make_agglom_cols(tax_glom_level = tax_glom_level, with_query_name = T)
  }

  if(glom_var == "n_unique_kmers"){
    taxonomy_annotate_df <- taxonomy_annotate_df %>%
      dplyr::select("genome_accession", "lineage", "query_name", "n_unique_kmers") %>%
      tidyr::separate(.data$lineage, into = c("domain", "phylum", "class", "order", "family", "genus", "species", "strain"), sep = ";", remove = F, fill = "right") %>%
      dplyr::group_by_at(dplyr::vars(dplyr::all_of(agglom_cols))) %>%
      dplyr::summarize(n_unique_kmers = sum(.data$n_unique_kmers)) %>%
      dplyr::ungroup() %>%
      tidyr::unite(col = "lineage", tidyselect::all_of(agglom_cols[-1]), sep = ";", remove = TRUE) %>%
      dplyr::select("lineage", "query_name", "n_unique_kmers")
  }
  # at the moment, this was easier to hard code with an if statement.
  # If I end up adding more glom vars, I'll figure out how to do this actually cleverly instead of copying and pasting the whole code chunk
  if(glom_var == "f_unique_to_query"){
    taxonomy_annotate_df <- taxonomy_annotate_df %>%
      dplyr::select("genome_accession", "lineage", "query_name", "f_unique_to_query") %>%
      tidyr::separate(.data$lineage, into = c("domain", "phylum", "class", "order", "family", "genus", "species", "strain"), sep = ";", remove = F, fill = "right") %>%
      dplyr::group_by_at(dplyr::vars(dplyr::all_of(agglom_cols))) %>%
      dplyr::summarize(f_unique_to_query = sum(.data$f_unique_to_query)) %>%
      dplyr::ungroup() %>%
      tidyr::unite(col = "lineage", tidyselect::all_of(agglom_cols[-1]), sep = ";", remove = TRUE) %>%
      dplyr::select("lineage", "query_name", "f_unique_to_query")
  }
  return(taxonomy_annotate_df)
}

# upset plot --------------------------------------------------------------

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

#' Visualize an upset plot of taxonomic lineage intersections between samples
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
    agglom_cols <- make_agglom_cols(tax_glom_level = tax_glom_level, with_query_name = F)

    # join together upset df with metadata in taxonomy_annotate_df
    taxonomy_annotate_df <-  taxonomy_annotate_df %>%
      dplyr::select("lineage") %>%
      tidyr::separate("lineage", into = agglom_cols, sep = ";", remove = F) %>%
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

# sankey diagram ----------------------------------------------------------

#' Visualize a sankey diagram from taxonomic lineages from one or many samples
#'
#' @description
#' `plot_taxonomy_annotate_sankey()` plots a sankey diagram from the output of sourmash taxonomy annotate.
#' The input data frame can contain one or many samples.
#' If there are many samples, abundances of each lineage are summarized and sample-level information is lost.
#' If the parameter `tax_glom_level` is specified, the plot will be summarized to that taxonomic rank (e.g. if "order" is specified, only domain, phylum, class, and order will be plotted).
#'
#' @param taxonomy_annotate_df Data frame containing outputs from sourmash taxonomy annotate. Can contain results from one or many runs of sourmash taxonomy annotate. If specified, agglomeration occurs across all queries.
#' @param tax_glom_level Optional character string specifying the taxonomic rank to agglomerate k-mer counts. Must be one of "domain", "phylum", "class", "order", "family", "genus", "species."
#' @param palette Optional character vector specifying a palette. Colors in the palette are recycled across taxonomic labels. If no palette is specified, RColorBrewer's Set2 is the default.
#'
#' @return A ggplot2 plot
#' @export
#'
#' @importFrom rlang .data
#'
#' @examples
#' \dontrun{
#' plot_taxonomy_annotate_sankey()
#' }
plot_taxonomy_annotate_sankey <- function(taxonomy_annotate_df, tax_glom_level = NULL, palette = NULL){
  if(!is.null(tax_glom_level)){
    agglom_cols <- make_agglom_cols(tax_glom_level = tax_glom_level, with_query_name = F)
  } else {
    agglom_cols <- c("domain", "phylum", "class", "order", "family", "genus", "species", "strain")
    # check if there are NAs in the strain column and emit a warning, as these will be dropped and the plot will look weird
    if(sum(is.na(taxonomy_annotate_df$strain)) > 0){
      stop("Some lineages are missing strain information. This will lead computation to fail for stat_parallel_sets_axes(). Use tax_glom_level = 'species' or a higher taxonomic rank to produce a plot")
    }
  }

  taxonomy_annotate_df <- taxonomy_annotate_df %>%
    tidyr::separate(.data$lineage, into = agglom_cols, sep = ";", remove = F, fill = "right", extra = "drop") %>%
    dplyr::group_by_at(dplyr::vars(tidyselect::all_of(agglom_cols))) %>%
    dplyr::summarize(sum_n_unique_kmers = sum(.data$n_unique_kmers))

  # format for ggforce parallel sets
  data <- ggforce::gather_set_data(taxonomy_annotate_df, 1:length(agglom_cols))

  # create a palette that recycles colors so each taxonomic label will be colorful
  if(is.null(palette)){
    # if the user doesn't supply a palette, use Set2
    palette <- c("#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3")
  }
  # otherwise ramp up from the user-defined palette
  palette <- grDevices::colorRampPalette(palette)(length(unique(data$y)))

  sankey_plt <- ggplot2::ggplot(data, ggplot2::aes(x = .data$x, id = .data$id, split = .data$y, value = .data$sum_n_unique_kmers)) +
    ggforce::geom_parallel_sets(alpha = 0.3, axis.width = 0.1) +
    ggforce::geom_parallel_sets_axes(axis.width = 0.2, ggplot2::aes(fill = .data$y)) +
    ggforce::geom_parallel_sets_labels(colour = 'black', angle = 360, size = 2, hjust = -0.25) +
    ggplot2::theme_classic() +
    ggplot2::theme(axis.line.y = ggplot2::element_blank(),
                   axis.text.y = ggplot2::element_blank(),
                   axis.ticks.y = ggplot2::element_blank(),
                   axis.ticks.x = ggplot2::element_blank(),
                   legend.position = "None") +
    ggplot2::labs(x = "tanomic rank") +
    # buffer the last axis so full names have space to print to viz
    ggplot2::scale_x_continuous(labels = c(agglom_cols, ""),
                                breaks = 1:(length(agglom_cols) + 1),
                                limits = c(.75, length(agglom_cols) + 1)) +
    ggplot2::scale_fill_manual(values = palette)

  return(sankey_plt)
}

# time series alluvial plot -----------------------------------------------

#' Visualize an allivual flow plot from taxonomic lineages from one or many samples
#'
#' @description
#' `plot_taxonomy_annotate_ts_alluvial()` creates an alluvial flow diagram from the output of sourmash taxonomy annotate for metagenomes that were sequenced in time series.
#'
#' @param taxonomy_annotate_df Data frame containing outputs from sourmash taxonomy annotate.
#' Can contain results from one or many runs of sourmash taxonomy annotate.
#' If specified, agglomeration occurs per query.
#' @param time_df A data frame. The first column should contain all of the values in the query_name column of taxonomy_annotate_df and should be called "query_name".
#' The second column should indicate the time that the sample was taken and should be named "time".
#' This plot requires that time be specified as numeric values (e.g. if samples were taken at 2 weeks, 1 month, and 2 months, the values should be specified as 0.5, 1, and 2);
#' this allows the function to appropriately sort and bound the x axis.
#' @param tax_glom_level Optional character string specifying the taxonomic rank to agglomerate the fraction of the metagenome that matched to a genome in the database (f_unique_to_query).
#' Must be one of "domain", "phylum", "class", "order", "family", "genus", "species."
#' @param fraction_threshold A number between 0-1. Defaults to 0.01.
#' The minimum fraction that a taxonomic lineage needs to occur in at least one time series sample for that lineage to have an alluvial ribbon in the final plot.
#' Lineages that occur below this threshold are grouped into an "other" category.
#'
#' @return A ggplot2 plot
#' @export
#'
#' @importFrom rlang .data
#'
#' @examples
#' \dontrun{
#' plot_taxonomy_annotate_ts_alluvial()
#' }
plot_taxonomy_annotate_ts_alluvial <- function(taxonomy_annotate_df, time_df, tax_glom_level = NULL, fraction_threshold = 0.01){
  # check formatting of time_df -- is the first column named query_name, and does it contain all of the query_names that are in the taxonomy_annotate_df?
  if(!all(colnames(time_df) == c("query_name", "time"))){
    stop("The column names of time_df must be query_name and time. Please update the column names using colnames(time_df) <- c('query_name', 'time') and re-run.")
  }

  if(!all(unique(taxonomy_annotate_df$query_name) %in% time_df$query_name)){
    stop("Not all query_name that are in taxonomy_annotate_df are in the query_name column of time_df. Please add the missing query_names and re-run")
  }

  # agglomerate to specified taxonomy level
  if(!is.null(tax_glom_level)){
    taxonomy_annotate_df <- tax_glom_taxonomy_annotate(taxonomy_annotate_df, tax_glom_level = tax_glom_level, glom_var = "f_unique_to_query")
  }

  # join the tax df with the time df
  taxonomy_annotate_df <- taxonomy_annotate_df %>%
    dplyr::left_join(time_df, by = c("query_name" = "query_name"))

  # determine which taxa to plot in their own alluvial ribbon.
  keep_taxa <- taxonomy_annotate_df %>%
    dplyr::filter(.data$f_unique_to_query > fraction_threshold)

  # determine agglom cols
  if(!is.null(tax_glom_level)){
    agglom_cols <- make_agglom_cols(tax_glom_level = tax_glom_level, with_query_name = F)
  } else {
    agglom_cols <- c("domain", "phylum", "class", "order", "family", "genus", "species", "strain")
    # check if there are NAs in the strain column and emit a warning, as these will be dropped and the plot will look weird
    if(sum(is.na(taxonomy_annotate_df$strain)) > 0){
      stop("Some lineages are missing strain information. This will lead to a very weird looking plot. Use tax_glom_level = 'species' or a higher taxonomic rank to produce a visualing rewarding plot")
    }
  }

  grp_by_vector <- c("query_name", "time", "tax_glom_col")
  alluvium_df <- taxonomy_annotate_df %>%
    tidyr::separate(.data$lineage, into = c("domain", "phylum", "class", "order", "family", "genus", "species", "strain"), sep = ";", remove = F, fill = "right") %>%
    # rename the column that will represent the alluvium ribbons to "tax_glom_col" so we can pass it to ggplot to as .data$
    dplyr::select("query_name", "time", "lineage", tax_glom_col = agglom_cols[length(agglom_cols)], "f_unique_to_query") %>%
    dplyr::mutate(tax_glom_col = ifelse(.data$lineage %in% keep_taxa$lineage, .data$tax_glom_col, "other")) %>%
    dplyr::select(-"lineage") %>%
    # creating the "other" designation creates a variable that is duplicated.
    # ggalluvium is not smart enough to sum over that variable itself a throws and error
    # the following lines of code won't change the value of f_unique_to_query for anything other than "other"
    dplyr::group_by_at(dplyr::vars(dplyr::all_of(grp_by_vector))) %>%
    dplyr::summarize(f_unique_to_query = sum(.data$f_unique_to_query))

  alluvial_plt <- ggplot2::ggplot(alluvium_df, ggplot2::aes(x = .data$time,
                                                  y = .data$f_unique_to_query,
                                                  alluvium = .data$tax_glom_col,
                                                  label    = .data$tax_glom_col,
                                                  fill     = .data$tax_glom_col)) +
    ggalluvial::geom_alluvium(colour = "black", alpha = .4, decreasing = FALSE) +
    ggplot2::labs(x = "Time", y = "Fraction of metagenome", colour = "", fill = tax_glom_level) +
    ggplot2::scale_y_continuous(expand = c(0, 0), limits = c(0, 1)) +
    ggplot2::theme_classic() +
    ggalluvial::stat_alluvium(geom = "text", size = 2, decreasing = FALSE, min.y = 0.005)

  return(alluvial_plt)
}

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
#' Inspired by phyloseq::tax_glom(), this method summarizes some numeric variables from genomes that have the same taxonomy at a user-specified taxonomy rank.
#' Agglomeration occurs within each sample, meaning the user-specified variable is only summed within each query_name.
#' This function returns a data frame with the columns `lineage`, `query_name`, and the glom_var column that was specified.
#' The accepted glom_vars that can be agglomerated are f_unique_to_query, f_unique_weighted, unique_intersect_bp, and n_unique_kmers.
#' Each of these variables deals with the "unique" fraction of the gather match, meaning there is no double counting between the query and the genome matched.
#' f_unique_weighted and n_unique_kmers are weighted by k-mer abundance while f_unique_to_query and unique_intersect_bp are not.
#' f_unique_weighted is similar to relative abundance (where "f" stands for fraction -- if 100% of the query had a match in the gather database, this value would sum to 1).
#' n_unique_kmers is the abundance-weighted number of unique hashes (k-mers) that overlapped between the query and the match in the database.
#' This number is calculated by dividing the unique_intersect_bp by the scaled value and multiplying this value by the average k-mer abundance.
#' Other variables (like f_orig_query) could sum > 1.
#' Only one variable is agglomerated at a time.
#'
#' @param taxonomy_annotate_df Data frame containing outputs from sourmash taxonomy annotate. Can contain results from one or many runs of sourmash taxonomy annotate. Agglomeration occurs within each query.
#' @param tax_glom_level Character. NULL by default, meaning no agglomeration is done. Valid options are "domain", "phylum", "class", "order", "family", "genus", and "species". When a valid option is supplied, k-mer counts are agglomerated to that level
#' @param glom_var Character. One of f_unique_to_query, f_unique_weighted, unique_intersect_bp, or n_unique_kmers.
#'
#' @return A data frame.
#' @export
#'
#' @importFrom rlang .data
#'
#' @details
#' Selecting which glom_var to use for downstream use cases can be difficult.
#' We most frequently use f_unique_weighted and n_unique_kmers as these both account for the number of times a k-mer occurs in a data set.
#' This is closer to counting the number of reads that would map against a reference genome than the other metrics.
#' When our downstream use case deals with relative abundance, f_unique_weighted is a good choice.
#' When the downstream use case needs count data, we use n_unique_kmers.
#' Because we divide by the scaled value to generate this number, the value will be much lower than read mapping.
#' However, doing it this way returns the actual number of k-mers sourmash counted.
#' This tends work better for assumptions made by downstream statistical tools (e.g. for differential abundance analysis, machine learning, etc.).
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

  if(!glom_var %in% c("f_unique_to_query", "f_unique_weighted", "unique_intersect_bp", "n_unique_kmers")){
    # these are the only ones that make sense to me to agglomerate across levels of taxonomy.
    stop("The variable you supplied is not a valid glom_var. Please choose from f_unique_to_query, f_unique_weighted, unique_intersect_bp, n_unique_kmers.")
  }

  taxonomy_annotate_df <- taxonomy_annotate_df %>%
    # temporarily rename the variable to "glom_var_tmp" so we can use the same code across all vars.
    dplyr::rename(glom_var_tmp = dplyr::all_of(glom_var)) %>%
    dplyr::select("genome_accession", "lineage", "query_name", "glom_var_tmp") %>%
    tidyr::separate(.data$lineage, into = c("domain", "phylum", "class", "order", "family", "genus", "species", "strain"), sep = ";", remove = F, fill = "right") %>%
    dplyr::group_by_at(dplyr::vars(dplyr::all_of(agglom_cols))) %>%
    dplyr::summarize(glom_var_tmp = sum(.data$glom_var_tmp)) %>%
    dplyr::ungroup() %>%
    tidyr::unite(col = "lineage", tidyselect::all_of(agglom_cols[-1]), sep = ";", remove = TRUE) %>%
    dplyr::select("lineage", "query_name", "glom_var_tmp")

  # rename the glom_var_tmp back to whatever the glom_var actually was
  colnames(taxonomy_annotate_df) <- c("lineage", "query_name", glom_var)

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
#' @param fill Optional argument specifying which level of taxonomy to fill the upset plot intersections with. Only levels above upset_inputs$tax_glom_level are valid. Uses the Set2 palette so cannot visualize more than 8 levels.
#' @param palette An optional character vector specifying the color palette to use. Ignored if fill is not set. Defaults to the colors in RColorBrewer Set2.
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
plot_taxonomy_annotate_upset <- function(upset_inputs, fill = NULL, palette = NULL){
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

    # create a palette if not user specified
    if(is.null(palette)){
      # if the user doesn't supply a palette, use Set2
      palette <- c("#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3")
    }

    # plot the upset plot
    plt <- ComplexUpset::upset(upset_df, intersect = unique(upset_inputs[[2]]$query_name), set_sizes = F,
                               base_annotations=list(
                                 '# lineages'=ComplexUpset::intersection_size(text=list(vjust=0.4, hjust=.05, angle=90),
                                                                              text_colors=c(on_background='black', on_bar='black'),
                                                                              mapping=ggplot2::aes(fill = .data$fill)) +
                                   ggplot2::scale_fill_manual(values = palette) +
                                   ggplot2::labs(fill = fill) +
                                   ggplot2::theme_classic() +
                                   ggplot2::theme(axis.text.x = ggplot2::element_blank(),
                                                  axis.ticks.x = ggplot2::element_blank(),
                                                  axis.title.x = ggplot2::element_blank()))
    )
    return(plt)
  }

  # plot the upset plot
  plt <- ComplexUpset::upset(upset_df, intersect = unique(upset_inputs[[2]]$query_name), set_sizes = F,
                             base_annotations=list(
                               '# lineages'=ComplexUpset::intersection_size(text=list(vjust=0.4, hjust=.05, angle=90),
                                                                            text_colors=c(on_background='black', on_bar='black')) +
                                 ggplot2::theme_classic() +
                                 ggplot2::theme(axis.text.x = ggplot2::element_blank(),
                                                axis.ticks.x = ggplot2::element_blank(),
                                                axis.title.x = ggplot2::element_blank())
  ))
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
    taxonomy_annotate_df <- tax_glom_taxonomy_annotate(taxonomy_annotate_df, tax_glom_level = tax_glom_level, glom_var = "f_unique_weighted")
  }

  # join the tax df with the time df
  taxonomy_annotate_df <- taxonomy_annotate_df %>%
    dplyr::left_join(time_df, by = c("query_name" = "query_name"))

  # determine which taxa to plot in their own alluvial ribbon.
  keep_taxa <- taxonomy_annotate_df %>%
    dplyr::filter(.data$f_unique_weighted > fraction_threshold)

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
    dplyr::select("query_name", "time", "lineage", tax_glom_col = agglom_cols[length(agglom_cols)], "f_unique_weighted") %>%
    dplyr::mutate(tax_glom_col = ifelse(.data$lineage %in% keep_taxa$lineage, .data$tax_glom_col, "other")) %>%
    dplyr::select(-"lineage") %>%
    # creating the "other" designation creates a variable that is duplicated.
    # ggalluvium is not smart enough to sum over that variable itself a throws and error
    # the following lines of code won't change the value of f_unique_to_query for anything other than "other"
    dplyr::group_by_at(dplyr::vars(dplyr::all_of(grp_by_vector))) %>%
    dplyr::summarize(f_unique_weighted = sum(.data$f_unique_weighted))

  # create a vector to use to reorder the legend so that "other" is always last
  if("other" %in% alluvium_df$tax_glom_col){
    order_vector <- alluvium_df %>%
      dplyr::filter(.data$tax_glom_col != "other")
    order_vector <- c(unique(order_vector$tax_glom_col), "other")
  } else {
    order_vector <- unique(alluvium_df$tax_glom_col)
  }

  alluvial_plt <- ggplot2::ggplot(alluvium_df, ggplot2::aes(x = .data$time,
                                                  y = .data$f_unique_weighted,
                                                  alluvium = .data$tax_glom_col,
                                                  label    = .data$tax_glom_col,
                                                  fill     = .data$tax_glom_col)) +
    ggalluvial::geom_alluvium(colour = "black", alpha = .4, decreasing = FALSE) +
    ggplot2::labs(x = "time", y = "abundance-weighted\nfraction of query", colour = "", fill = tax_glom_level) +
    ggplot2::scale_y_continuous(expand = c(0, 0), limits = c(0, 1)) +
    ggplot2::theme_classic() +
    ggalluvial::stat_alluvium(geom = "text", size = 2, decreasing = FALSE, min.y = 0.005) +
    ggplot2::scale_fill_discrete(breaks = order_vector)

  return(alluvial_plt)
}

# predict number of strains for each species detected ---------------------

#' Detect whether multiple strains of the same species are present in a sample.
#'
#' @description
#' `from_taxonomy_annotate_to_multi_strains()` uses the output of sourmash taxonomy annotate to detect whether multiple strains of the same species are present in a sample (e.g. a metagenome).
#' The function uses the `f_match` metric produced by sourmash gather to predict if multiple strains are present.
#' `f_match` is the fraction of a matched genome that is contained within the query sample.
#' This function sums over all `f_match` values for all matched genomes of a given species and detects when we see more genomic segments than we would expect to see if only one strain were truly present.
#' Most frequently, if multiple genomes form one species are identified by sourmash gather, the single genome in a query sample is different from genomes in databases so multiple genomes that each cover a distinct portion of the genome that is truly present are returned as matches.
#' However, for each species, the fraction of a single genome that matched against the query decreases with each successive match.
#' This is consistent with the idea of pangenomes -- the single true genome in our query matches many different parts of genomes in the database.
#' When this happens, the sum fraction matched (`f_match`) by all genomes within a species within a single query should not exceed ~1.
#' When we run sourmash gather on single genomes, the real value we observe for the sum f_match ranged between:
#' 0.046 (for genomes that only had ~genus-level relatives in the database) and 1.04 (for genomes that had many close matches in the database).
#' We feel reasonably confident that it is a good starting place for strain-level analysis (see details below for caveats) but additional validated methods should be used to confirm findings.
#'
#' @details
#' In this context, we refer to _strain_ as any sub-species level variation.
#' This function uses genomes as the unit of strain -- each genome is considered a different strain.
#' Sourmash gather compares a query (e.g. a metagenome) against a database (such as GenBank microbial genomes) and provides the minimum set of genomes that cover all of the k-mers in the query that are in the database.
#' At least two things could be happening when sourmash gather returns two genomes of the same species (e.g. different strains) as a match to the same metagenome sample:
#' 1. Both strains may be present in the metagenome
#' 2. Only one strain may be truly present in the metagenome, but that strain is not contained within our current reference database.
#' Instead, pieces of that strain's genome are in other genomes in the database.
#' The genome in the database that contains the largest overlap with the genome in the metagenome is returned first as the best match.
#' Then, other genomes in the database are returned that match other portions of the metagenome strain's genome that wasn't contained in the best match.
#' In reality, some combination of these two things probably happens.
#' Keep in mind:
#' 1. We can't detect strain variation if there is only one genome for a given species in the database using sourmash gather/taxonomy alone (variant calling tools or tools that look at variation in assembly graphs would be more successful for this use case).
#' 2. The sourmash gather/taxonomy results alone should not be used to conclusively detect strain variation, but it is a good place to start to figure out where to dig in deeper.
#'
#' @param taxonomy_annotate_df Data frame containing outputs from sourmash taxonomy annotate.
#' Can contain results from one or many runs of sourmash taxonomy annotate.
#' @param plot_threshold f_match threshold for plotting a genome match.
#' This threshold is for plotting only.
#'
#' @return A named list.
#' The first object in the list `candidate_species_with_multiple_strains` summarizes the query sample and species names that may have an f_match >= 1.1.
#' The second object `plt` is a ggplot object that summarizes the query samples and species that may contain multiple strains.
#' The third object `plt_data` contains the data that is used to produce the ggplot object.
#' @export
#'
#' @importFrom rlang .data
#'
#' @examples
#' \dontrun{
#' from_taxonomy_annotate_to_multi_strains()
#' }
from_taxonomy_annotate_to_multi_strains <- function(taxonomy_annotate_df, plot_threshold = 0.02){
  # check if query name is all NAs, fill it in with query_filename
  if(all(is.na(taxonomy_annotate_df$query_name))){
    taxonomy_annotate_df <- taxonomy_annotate_df %>%
      dplyr::mutate(query_name = ifelse(is.na(.data$query_name), basename(.data$query_filename), .data$query_name))
  }

  # for each query in taxonomy_annotate_df, count how many genomes are observed per species
  more_than_one_genome_observed_for_species <- taxonomy_annotate_df %>%
    dplyr::group_by(.data$query_name, .data$species) %>%
    dplyr::tally() %>%
    dplyr::filter(.data$n > 1)

  # sum the f_match within each sample and species
  f_match <- taxonomy_annotate_df %>%
    dplyr::filter(.data$species %in% more_than_one_genome_observed_for_species$species) %>% # filter to species with more than one genome observed
    dplyr::group_by(.data$query_name, .data$species) %>%
    dplyr::summarise(species_f_match = sum(.data$f_match)) %>%
    dplyr::arrange(dplyr::desc(.data$species_f_match))

  # filter to species that summed to an f_match >= 1.1
  f_match_filtered <- f_match %>%
    dplyr::filter(.data$species_f_match >= 1.1) %>%
    # create a vector to filter with so we only get combinations that passed the filter
    dplyr::mutate(query_name_species = paste0(.data$query_name, "-", .data$species))

  # check if there were any matches. If no matches, print a helpful message, return an empty list of the same structure as what is returned if there are results, and exit
  if(nrow(f_match_filtered) == 0){
    return(list(candidate_species_with_multiple_strains = f_match_filtered,
                plt = NULL,
                plt_data = NULL))
    stop("There were no species with potential multiple strains in the supplied metagenomes given the database used during sourmash gather and sourmash taxonomy.")
  }

  # filter to the species that meet the criteria and to show in the plot
  plt_df <- taxonomy_annotate_df %>%
    dplyr::mutate(query_name_species = paste0(.data$query_name, "-", .data$species)) %>%
    dplyr::filter(.data$query_name_species %in% f_match_filtered$query_name_species) %>%
    # for plotting only, filter to genomes that have an f_match of plot_threshold or greater
    dplyr::filter(.data$f_match >= plot_threshold) %>%
    dplyr::mutate(species = paste0("italic('", .data$species, "')"))

  # produce a plot
  plt <- ggplot2::ggplot(plt_df, ggplot2::aes(x = stats::reorder(.data$genome_accession, -.data$f_match),
                                              y = .data$average_abund,
                                              label = round(.data$f_match, digits = 2))) +
    ggplot2::geom_point(ggplot2::aes(size = .data$f_match)) +
    ggplot2::coord_flip() +
    ggplot2::facet_wrap(~.data$query_name + .data$species, scales = "free",
                        labeller = ggplot2::label_parsed) +
    ggrepel::geom_text_repel(size = 2, color = "grey") +
    ggplot2::theme_classic() +
    ggplot2::theme(strip.background = ggplot2::element_blank()) +
    ggplot2::labs(y = "average abundance of k-mers in genome",
                  x = "genome",
                  size = "fraction of genome")

  # return a list with plot and summaries
  return(list(candidate_species_with_multiple_strains = f_match_filtered,
              plt = plt,
              plt_data = plt_df))
}

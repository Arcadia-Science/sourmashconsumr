
# fraction of sample identified bar char ----------------------------------

#' Visualize the classified and unclassified fraction of many samples using sourmash gather results.
#'
#' @description
#' `plot_gather_classified()` produces a ggplot2 plot of the classified and unclassified fraction of many samples from a data frame of sourmash gather results.
#'
#' @param gather_df A data frame of multiple sourmash gather results created by `read_gather()`.
#'
#' @return A ggplot2 plot.
#' @export
#'
#' @importFrom rlang .data
#'
#' @examples
#' \dontrun{
#' plot_gather_classified(gather_df)
#' }
plot_gather_classified <- function(gather_df) {
  # generate a dataframe that contains information about the percent of the sample that was unclassified.
  # then, bind this new dataframe to the sourmash_taxonomy dataframe within the plot
  # that way, this information is not carried on beyond this plot.
  gather_df <- read_gather(Sys.glob("tests/testthat/*gather*.csv")) %>%
    dplyr::mutate(database = basename(filename))

  gather_df_unclassified <- gather_df %>%
    dplyr::group_by(.data$query_name) %>%
    dplyr::summarize(sum_f_unique_weighted = sum(.data$f_unique_weighted)) %>%
    dplyr::mutate(f_unique_weighted = 1 - sum_f_unique_weighted,
                  database = "unclassified") %>%
    dplyr::select(-"sum_f_unique_weighted")

  gather_df <- gather_df %>%
    dplyr::bind_rows(gather_df_unclassified)

  ggplot2::ggplot(gather_df,
                  ggplot2::aes(x = query_name, y = f_unique_weighted, fill = database)) +
    ggplot2::geom_col() +
    ggplot2::theme_classic() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90)) +
    ggplot2::scale_fill_brewer(palette = "Paired") +
    ggplot2::labs(x = "query name", y = "abundance-weighted\nfraction of query")
}

# upset plots -------------------------------------------------------------

#' Convert a data frame containing sourmash gather results into a upset plot compatible data frame
#'
#' @description
#' `from_gather_to_upset_df()` converts a data frame like that produced by `read_gather()` into a named list of vectors.
#' Each vector in the list is named by the query_name column in the input data frame and contains the genome accessions found in that sample.
#' The `query_name` column is used for sample names.
#' The returned data frame will have the genome accessions as row names and the query_names column names and can be plotted with `UpSetR::upset()` or `ComplexUpset::upset()`.
#' If plotting with Complex Upset, additional metadata can be added to the data frame (and therefore the plot) by joining on the values of the genome_accession rownames.
#'
#' @param gather_df A data frame of multiple sourmash gather results created by `read_gather()`.
#'
#' @return An upset plot compliant data frame
#' @export
#'
#' @importFrom rlang .data
#'
#' @examples
#' \dontrun{
#' from_gather_to_upset_df(gather_df)
#' }
from_gather_to_upset_df <- function(gather_df){
  # turn dataframe into list.
  # Each index in the list is named after the signature it represents.
  # Each index contains a vector of genome accessions that were in that signature
  . <- NULL # need to set . as global var because of whack line.
  # I don't actually like this approach, but I think it's ok for now until I can find a workaround for said line.
  upset_list <- gather_df %>%
    dplyr::group_by(.data$query_name) %>%
    {stats::setNames(dplyr::group_split(.), dplyr::group_keys(.)[[1]])} # absolutely whack line of code that allows us to set the list names by the group_by variable

  # extract just the minhashes
  upset_list <- lapply(upset_list, function(x) {x$genome_accession})

  # convert list into upset-compliant data.frame
  upset_df <- from_list_to_upset_df(upset_list)

  return(upset_df)
}

#' Visualize the intersection of genome accessions in sourmash gather results from many samples
#'
#' @description
#' `plot_gather_upset` uses `ComplexUpset::upset()` to plot the intersection of genome accessions in sourmash gather results from many samples.
#' The plot can optionally be colored by the database from which the genome accession originated.
#' While based on ggplot2, `ComplexUpset::upset()` has a difficult syntax to parameterize.
#' To produce an alternative visualization, it may be easiest to run `sourmashconsumr::plot_gather_upset` to retrieve the function source code and alter it.
#'
#' @param upset_df An upset plot compliant data frame like that produced by from_gather_to_upset_df.
#' @param color_by_database Boolean indicating whether to fill the upset plot instersection bar plot by the database the results came from. FALSE by default.
#' @param gather_df Gather results passed to from_gather_to_upset_df to create the input data frame. NULL unless color_by_database is set to TRUE.
#'
#' @return A ComplexUpset plot.
#' @export
#'
#' @importFrom rlang .data
#'
#' @examples
#' \dontrun{
#' plot_gather_upset(upset_df)
#' }
plot_gather_upset <- function(upset_df, color_by_database = FALSE, gather_df = NULL){
  if(color_by_database == F){
  upset_plt <- ComplexUpset::upset(upset_df, intersect = names(upset_df), set_sizes = F,
                                   base_annotations=list(
                                     '# scaled k-mers' = ComplexUpset::intersection_size(text=list(vjust=0.4, hjust=.05, angle=90),
                                                                                       text_colors=c(on_background='black', on_bar='black'),
                                                                                       mapping=ggplot2::aes(fill='bars_color')) +
                                       ggplot2::scale_fill_manual(values=c('bars_color'='lightgrey'), guide='none'))
  )
  }
  if(color_by_database == T){
    # make sure gather_df is not NULL
    if(is.null(gather_df)){
      stop("gather_df must be defined to color_by_database. Please use the gather_df you used as input to from_gather_to_upset_df().")
    }

    db_df <- gather_df %>%
      dplyr::select('genome_accession', 'database') %>%
      dplyr::distinct()

    # check and make sure that there aren't repeats (e.g. the same accession in multiple databases)
    # this will fail without me putting a stop in (bc there there can't be duplicate row names)
    # but with the stop it will fail informatively so the user can fix the problem themselves
    check_db_df <- db_df %>%
      dplyr::group_by(.data$genome_accession) %>%
      dplyr::tally()
    if(any(check_db_df$n > 1)){
      stop("At least one genome accession was in multiple databases for the supplied gather results. Either use color_by_database == FALSE or limit which samples you plot.")
    }

    upset_df <- upset_df %>%
      tibble::rownames_to_column("genome_accession") %>%
      dplyr::left_join(db_df, by = "genome_accession") %>%
      tibble::column_to_rownames("genome_accession")

    ComplexUpset::upset(upset_df, intersect = names(upset_df)[1:(ncol(upset_df)-1)], set_sizes = F,
                        base_annotations = list(
                          '# lineages' = ComplexUpset::intersection_size(text=list(vjust=0.4, hjust=.05, angle=90),
                                                                         text_colors=c(on_background='black', on_bar='black'),
                                                                         mapping=ggplot2::aes(fill = .data$database)) +
                            ggplot2::scale_fill_brewer(palette = "Set2"))
    )
  }
  return(upset_plt)
}

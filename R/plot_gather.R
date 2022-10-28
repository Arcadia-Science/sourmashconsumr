
# fraction of sample identified bar char ----------------------------------

plot_gather_classified <- function(gather_df) {
  # generate a dataframe that contains information about the percent of the sample that was unclassified.
  # then, bind this new dataframe to the sourmash_taxonomy dataframe within the plot
  # that way, this information is not carried on beyond this plot.
  gather_df <- read_gather(Sys.glob("tests/testthat/*gather*.csv")) %>%
    dplyr::mutate(database = basename(filename))
  gather_df_unclassified <- gather_df %>%
    dplyr::group_by(query_name) %>%
    dplyr::summarize(sum_f_unique_weighted = sum(f_unique_weighted)) %>%
    dplyr::mutate(f_unique_weighted = 1 - sum_f_unique_weighted,
                  database = "unclassified") %>%
    dplyr::select(-sum_f_unique_weighted)

  ggplot2::ggplot(gather_df %>%
           dplyr::bind_rows(gather_df_unclassified),
         ggplot2::aes(x = query_name, y = f_unique_weighted, fill = database)) +
    ggplot2::geom_col() +
    ggplot2::theme_classic() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90)) +
    ggplot2::scale_fill_brewer(palette = "Paired")
}

# upset plots -------------------------------------------------------------

#' Title
#'
#' @param gather_df
#'
#' @return An upset plot compliant data frame
#' @export
#'
#' @importFrom rlang .data
#'
#' @examples
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

plot_gather_upset <- function(upset_df, color_by_database = T, gather_df){
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

    upset_df <- upset_df %>%
      tibble::rownames_to_column("genome_accession") %>%
      dplyr::left_join(db_df, by = "genome_accession") %>%
      tibble::column_to_rownames("genome_accession")

    ComplexUpset::upset(upset_df, intersect = names(upset_df), set_sizes = F,
                        base_annotations = list(
                          '# lineages' = ComplexUpset::intersection_size(text=list(vjust=0.4, hjust=.05, angle=90),
                                                                         text_colors=c(on_background='black', on_bar='black'),
                                                                         mapping=ggplot2::aes(fill = .data$database)) +
                            ggplot2::scale_fill_brewer(palette = "Set2") +
                            ggplot2::labs(fill = fill))
    )
  }
  return(upset_plt)
}

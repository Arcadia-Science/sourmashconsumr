#' Internal function to check if the sample name is a column or a rowname and move it to a rowname if it is a column.
#'
#' @param compare_df
#'
#' @return A data frame.
#'
#' @examples
#' check_compare_df_sample_col_and_move_to_rowname(compare_df)
check_compare_df_sample_col_and_move_to_rowname <- function(compare_df){
  # check if compare data frame was read in with sample names as a column
  if("sample" %in% colnames(compare_df)){
    # move sample column to rownames
    compare_df <- compare_df %>%
      tibble::column_to_rownames("sample")
  }
  return(compare_df)
}


#' Perform an MDS analysis from a sourmash compare output and produce a tidy data.frame
#'
#' @param compare_df A data frame produced by read_compare_csv() representing a (dis)similarity matrix output by sourmash compare.
#'
#' @return A data frame.
#' @export
#'
#' @importFrom rlang .data
#'
#' @examples
#' make_compare_mds(compare_df)
make_compare_mds <- function(compare_df){
  # check if compare tibble was read in with sample names as a column
  compare_df <- check_compare_df_sample_col_and_move_to_rowname(compare_df)

  compare_mds <- compare_df %>%
    as.matrix() %>%
    stats::dist() %>%
    stats::cmdscale() %>%
    as.data.frame() %>%
    dplyr::rename(MDS1 = .data$V1, MDS2 = .data$V2) %>%
    tibble::rownames_to_column("sample")

  return(compare_mds)
}

#' Plot an MDS data.frame produced from the output of sourmash compare
#'
#' @param compare_mds A data.frame produced using make_compare_mds()
#' @param label Boolean controlling whether sample labels are added to the plot. The default, TRUE, plots labels.
#'
#' @return A ggplot2 plot.
#' @export
#'
#' @importFrom rlang .data
#'
#' @examples
#' plot_compare_mds(compare_mds)
plot_compare_mds <- function(compare_mds, label = TRUE){
  mds_plt <- ggplot2::ggplot(compare_mds, ggplot2::aes(x = .data$MDS1, y = .data$MDS2, label = .data$sample)) +
    ggplot2::geom_point() +
    ggplot2::theme_classic()

  if(label == TRUE){
    mds_plt <- mds_plt +
      ggrepel::geom_label_repel()
  }

  return(mds_plt)
}

#' Plot a heatmap from a tibble of sourmash compare results.
#'
#' @param compare_df A tibble or data.frame produced by read_compare_csv() representing a (dis)similarity matrix output by sourmash compare.
#' @param seed Random seed for function. Default is 42. Ensures repeatable plot outputs.
#' @param ... Arguments passed to base::heatmap().
#'
#' @return A plot of sourmash compare results.
#' @export
#'
#' @examples
#' plot_compare_heatmap()
plot_compare_heatmap <- function(compare_df, seed = 42, ...){
  # set seed so the same visualization is produced each time the code is run
  set.seed(seed = seed)

  # check if compare tibble was read in with sample names as a column
  compare_df <- check_compare_df_sample_col_and_move_to_rowname(compare_df)

  compare_mat <- as.matrix(compare_df)

  # use hclust() and dist() to calculate clusters
  hc_rows <- stats::hclust(stats::dist(compare_mat), method = "single") # method single matches sourmash plot
  hc_cols <- stats::hclust(stats::dist(t(compare_mat)), method = "single")

  # plot the heatmap
  stats::heatmap(compare_mat,
                 Colv = stats::as.dendrogram(hc_cols),
                 Rowv = stats::as.dendrogram(hc_rows),
                 scale = 'none', ...)
}

make_compare_mds <- function(){

}

plot_compare_mds <- function(){}

#' Plot a heatmap from a tibble of sourmash compare results.
#'
#' @param compare_df A tibble containing results produced by sourmash.
#' @param ... Arguments passed to base::heatmap().
#'
#' @return A plot of sourmash compare results.
#' @export
#'
#' @examples
#' plot_compare_heatmap()
plot_compare_heatmap <- function(compare_df, ...){
  # check if compare tibble was read in with sample names as a column
  if("sample" %in% colnames(compare_df)){
    # guard against missing pkg installations
    if (!requireNamespace("tibble", quietly = TRUE)) {
      stop(
        "Package \"tibble\" must be installed to make sample names into row names. Please install \"tibble\".",
        call. = FALSE
      )
    }
    # move sample column to rownames
    compare_df <- compare_df %>%
      tibble::column_to_rownames("sample")
  }

  compare_mat <- as.matrix(compare_df)

  # use hclust() and dist() to calculate clusters
  hc_rows <- hclust(dist(compare_mat), method = "single") # method single matches sourmash plot
  hc_cols <- hclust(dist(t(compare_mat)), method = "single")

  # plot the heatmap
  heatmap(compare_mat,
          Colv = as.dendrogram(hc_cols),
          Rowv = as.dendrogram(hc_rows),
          scale='none', ...)
}


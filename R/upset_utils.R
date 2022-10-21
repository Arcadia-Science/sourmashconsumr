#' Convert a list of vectors into an upset plot-compatible data.frame while conserving value names as row.names
#'
#' @description
#' `from_list_to_upset_df()` copies the functionality of UpSetR::fromList() and converts a list of named vectors to an upset plot compatible data frame.
#' It modifies the original function to return a data frame with row names that record the identity of the items in each vector.
#' This allows additional metadata to be joined to the resultant data frame using the values of the items in the list.
#'
#' @param upset_list A named list of vectors. Each name should be one set and each vector to should contain the values to be intersected.
#'
#' @return A data.frame.
#' @export
#'
#' @examples
#' sample1 <- c("a", "b", "c")
#' sample2 <- c("b", "c", "d")
#' upset_list <- list(sample1 = sample1, sample2 = sample2)
#' from_list_to_upset_df(upset_list)
from_list_to_upset_df <- function (upset_list) {
  # source: https://github.com/hms-dbmi/UpSetR/issues/85#issuecomment-327900647
  # Same as original UpSetR::fromList()...
  # function creates an upset data frame where each column represents a sample and each row represents an observation.
  # Rows are values 0/1 for presence/absence of the observation within the sample.
  # the behavior of this function differs from UpSetR::fromList() in that this one preserves the names of the observations as row names in the final DF.
  # this is helpful for ComplexUpset plots so that the intersection bar charts can be colored by metadata associated with the observation.
  # Note that UpSetR::fromList() is under an MIT license (checked October 2022).
  elements <- unique(unlist(upset_list))
  data <- unlist(lapply(upset_list, function(x) {
    x <- as.vector(match(elements, x))
  }))
  data[is.na(data)] <- as.integer(0)
  data[data != 0] <- as.integer(1)
  data <- data.frame(matrix(data, ncol = length(upset_list), byrow = F))
  data <- data[which(rowSums(data) != 0), ]
  names(data) <- names(upset_list)
  # ... Except now it conserves your original value names!
  row.names(data) <- elements
  return(data)
}

#' Transform an upset compliant data frame into a data frame of intersection names and lineages
#'
#' @description
#' `from_upset_df_to_intersections()` transforms an upset compliant data frame into a data frame that details the intersections shared by the samples in the upset data frame.
#' It produces a data frame with a column named `intersection`.
#' If the intersection is shared between a set of samples, their column names occur in the character string in the `intersection` column.
#' The other sample column names appear as 0s.
#' The identity of the intersection is recorded in the row names.
#'
#' @param upset_df An upset compliant data frame.
#'
#' @return A data frame.
#' @export
#'
#' @examples
#' sample1 <- c("a", "b", "c")
#' sample2 <- c("b", "c", "d")
#' upset_list <- list(sample1 = sample1, sample2 = sample2)
#' upset_df <- from_list_to_upset_df(upset_list)
#' from_upset_df_to_intersections(upset_df)
from_upset_df_to_intersections <- function(upset_df){
  # if the upset_df came out of from_taxonomy_annotate_to_upset_inputs() or from_signatures_to_upset_df(), the upset df will only contain 0/1 columns from samples.
  # the lineage information will be in the row name
  # replace 1s in the upset_df with names of column.
  # this determines the number of columns in the input data frame
  intersection_df <- data.frame(lapply(seq_along(upset_df),
                                       function(x) ifelse(upset_df[[x]] == 1,
                                                          names(upset_df)[x],
                                                          upset_df[[x]])))
  # reset dropped rownames
  rownames(intersection_df) <- rownames(upset_df)
  # create a new column in the data frame named "intersection" that records the the groups of samples that share a given lineage
  intersection_df <- intersection_df %>%
    tidyr::unite(col = "intersection", sep = "_")

  return(intersection_df)
}

#' Summarize the size of intersections between samples in an upset compliant data frame
#'
#' @description
#' `from_upset_df_to_intersection_summary()` reports the same information as upset plots but in a data frame.
#' It reports the intersection identity and the number of elements shared by the intersections.
#'
#' @param upset_df An upset compliant data frame.
#'
#' @return A data frame.
#' @export
#'
#' @importFrom rlang .data
#'
#' @examples
#' sample1 <- c("a", "b", "c")
#' sample2 <- c("b", "c", "d")
#' upset_list <- list(sample1 = sample1, sample2 = sample2)
#' upset_df <- from_list_to_upset_df(upset_list)
#' from_upset_df_to_intersection_summary(upset_df)
from_upset_df_to_intersection_summary <- function(upset_df){
  intersection_df <- from_upset_df_to_intersections(upset_df)
  # create a data frame showing intersections and the number of shared lineages
  intersection_df_summarized <- intersection_df %>%
    dplyr::group_by(.data$intersection) %>%
    dplyr::summarise(n = dplyr::n()) %>%
    dplyr::arrange(dplyr::desc(.data$n))

  return(intersection_df_summarized)
}

#' Extract the lineages that underlie the intersections in an upset plot
#'
#' @description
#' `from_upset_df_to_intersection_members()` produces a named list where each item in the list a vector of the items that make up an intersection.
#' The names of each item in the list are concatenations of sample names (built from `colnames(upset_df)`).
#' If the sample name occurs in the name, the members in the vectors were observed in that sample name.
#' Zeroes (`0`) are used as place holders for samples that don't have a given member.
#'
#' @param upset_df An upset compliant data frame.
#'
#' @return A list.
#' @export
#'
#' @examples
#' sample1 <- c("a", "b", "c")
#' sample2 <- c("b", "c", "d")
#' upset_list <- list(sample1 = sample1, sample2 = sample2)
#' upset_df <- from_list_to_upset_df(upset_list)
#' lst <- from_upset_df_to_intersection_members(upset_df)
#' # to pull out values from the list, index into the list using any of the intersection names.
#' lst[["sample1_sample2"]] # identity of values shared by sample1 and sample2
#' lst["sample1_0"] # identity of values that occurred in only sample1
from_upset_df_to_intersection_members <- function(upset_df){
  intersection_df <- from_upset_df_to_intersections(upset_df)

  intersection_members_list <- intersection_df %>%
    tibble::rownames_to_column("lineage") %>%
    dplyr::group_by(.data$intersection) %>%
    dplyr::summarise(list = list(.data$lineage)) %>%
    dplyr::mutate(list = stats::setNames(list, .data$intersection)) %>%
    dplyr::pull(list)

  return(intersection_members_list)
}

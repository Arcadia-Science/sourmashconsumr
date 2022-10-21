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

summarize_upset_df_intersections <- function(upset_df){}

list_upset_df_intersections <- function(upset_df){}

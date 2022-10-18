#' Convert a list of vectors into an upset plot-compatible data.frame while conserving value names as row.names
#'
#' @param list
#'
#' @return
#' @export
#'
#' @examples
from_list_to_upset_df <- function (upset_list) {
  # source: https://github.com/hms-dbmi/UpSetR/issues/85#issuecomment-327900647
  # Same as original UpSetR::fromList()...
  # function creates an upset data frame where each column represents a sample and each row represents an observation.
  # Rows are values 0/1 for presence/absence of the observation within the sample.
  # the behavior of this function differs from UpSetR::fromList() in that this one preserves the names of the observations as row names in the final DF.
  # this is helpful for ComplexUpset plots so that the intersection bar charts can be colored by metadata associated with the observation.
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

#' Check uniformity in sketch parameters
#'
#' @param signatures A data.frame representing many signatures.
#'
#' @return NULL. Stops execution if parameters are not uniform.
#'
#' @examples
#' check_uniform_parameters_in_signatures_df(signatures)
check_uniform_parameters_in_signatures_df <- function(signatures){
  # check that there is only one k-mer size, scaled value, num, hash_function, seed, molecule type in the combined data.frame.
  # It doesn't make sense to visualize across these values as they represent fundamentally different things that are not inter-operable.

  check_ksize <- length(unique(signatures$ksize))
  check_molecule <- length(unique(signatures$molecule))
  check_hash_function <- length(unique(signatures$hash_function))
  check_seed <- length(unique(signatures$seed))
  check_scaled <- length(unique(signatures$seed))
  check_num <- length(unique(signatures$seed))

  if(check_ksize != 1){
    stop("The signatures represented in this data frame have different k-mer sizes. Please filter to one ksize before continuing.")
  }

  if(check_molecule != 1){
    stop("The signatures represented in this data frame have different molecule types. Please filter to one molecule type before continuing.")
  }

  if(check_hash_function != 1){
    stop("The signatures represented in this data frame were sketched with different hash functions. Please filter to one hash function before continuing.")
  }

  if(check_seed != 1){
    stop("The signatures represented in this data frame were sketched with different seeds. Please filter to one seed before continuing.")
  }

  if(check_scaled != 1 & check_num != 1){
    # only scaled or num can be set, so check at the same time.
    stop("The signatures represented in this data frame were sketched with scaled or num values. Please filter to one scaled or num value before continuing.")
  }
}

#' Convert a data.frame containing information from many sourmash signatures into a upset plot-compatible data.frame
#'
#' @param signatures A data.frame of multiple sourmash signatures created by combining many signatures read in by read_signature()
#'
#' @return
#' @export
#'
#' @examples
from_signatures_to_upset_df <- function(signatures){
  # stop if not all parameters to build the sketches are uniform, otherwise the intersections in the upset plot won't make sense
  check_uniform_parameters_in_signatures_df(signatures)

  # turn dataframe into list.
  # Each index in the list is named after the signature it represents.
  # Each index contains a vector of hashes that were in that signature
  upset_list <- signatures %>%
    dplyr::mutate(name = basename(filename)) %>% # build names based on file name
    dplyr::group_by(name) %>%
    {setNames(dplyr::group_split(.), dplyr::group_keys(.)[[1]])} # absolutely whack line of code that allows us to set the list names by the group_by variable

  # extract just the minhashes
  upset_list <- lapply(upset_list, function(x) {x$mins})

  # convert list into upset-compliant data.frame
  upset_df <- from_list_to_upset_df(upset_list)

  return(upset_df)
}


#' Visualize the intersection of hashes in sourmash signatures from many samples
#'
#' @param upset_df An upset-compliant data.frame.
#'
#' @return
#' @export
#'
#' @examples
plot_signatures_upset <- function(upset_df){
  upset_plt <- ComplexUpset::upset(upset_df, intersect = names(upset_df), set_sizes = F,
                                   base_annotations=list(
                                     '# scaled k-mers'=ComplexUpset::intersection_size(text=list(vjust=0.4, hjust=.05, angle=90),
                                                                                       text_colors=c(on_background='black', on_bar='black'),
                                                                                       mapping=ggplot2::aes(fill='bars_color')) +
                                       ggplot2::scale_fill_manual(values=c('bars_color'='lightgrey'), guide='none'))
  )
  return(upset_plt)
}

#' Create a tidy data.frame representing per-sample rarefaction information
#'
#' @param signatures A data.frame of multiple sourmash signatures created by combining many signatures read in by read_signature()
#'
#' @return
#' @export
#'
#' @examples
from_signatures_to_rarefaction_df <- function(signatures){
  # check that the signatures dataframe has a column for abundances
  if("abundances" %in% colnames(signatures) == FALSE){
    stop("The signatures data.frame does not have an 'abundances' columns.
         Abundances are required to calculate rarefaction cureves.
         If you have abundances in your signatures data.frame, make sure the column is named 'abundances'.")
  }

  # format the signatures into a rarecurve-compliant data.frame
  signatures <- signatures %>%
    dplyr::mutate(name = basename(filename)) %>% # set names as basename of filename
    dplyr::select(name, mins, abundances) %>% # select columns with relevant data
    tidyr::pivot_wider(id_cols = name, values_from = "abundances", names_from = "mins") %>% # transform to wide format
    tibble::column_to_rownames("name") %>%
    replace(is.na(.), 0) # back fill NAs with 0s

  # remove hashes (columns) that are only observed once and in one sample
  # since these rarefaction curves should only be built on reads, these hashes are probably errors
  signatures <- signatures %>%
    dplyr::select(which(!colSums(signatures, na.rm=TRUE) %in% 1))

  # return as a tidy data frame with per-sample information
  # this data frame will include the information needed to plot the rarefaction curves
  rarecurve_df <- vegan::rarecurve(signatures, step = 2, tidy = T)

  # make the column names more meaningful
  colnames(rarecurve_df) <- c("name", "num_kmers_sampled", "num_kmers_observed")

  return(rarecurve_df)
}

#' Plot a rarefaction curve to determine sequence saturation
#'
#' @param rarefaction_df
#' @param fraction_of_points_to_plot
#'
#' @return
#' @export
#'
#' @examples
plot_signatures_rarefaction <- function(rarefaction_df, fraction_of_points_to_plot = 500){
  # we don't need to plot all of the points on the line to get an accurate picture of sequencing saturation
  # This filters to 1/fraction_of_points_to_plot at the correct intervals to maintain curve shape
  rarefaction_df_filtered <- rarefaction_df %>%
    dplyr::filter(num_kmers_sampled %% fraction_of_points_to_plot == 1)

  plt <- ggplot2::ggplot(rarefaction_df_filtered,
                ggplot2::aes(x = num_kmers_sampled, y = num_kmers_observed)) +
    ggplot2::geom_point() +
    ggplot2::labs(x = "number of sampled k-mers", y = "number of distinct k-mers observed")

  return(plt)
}

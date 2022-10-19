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

#' Check uniformity in sketch parameters
#'
#' @param signatures A data.frame representing many signatures.
#'
#' @return NULL. Stops execution if parameters are not uniform.
#'
#' @examples
#' check_uniform_parameters_in_signatures_df(signatures)
check_uniform_parameters_in_signatures_df <- function(signatures_df){
  # check that there is only one k-mer size, scaled value, num, hash_function, seed, molecule type in the combined data.frame.
  # It doesn't make sense to visualize across these values as they represent fundamentally different things that are not inter-operable.

  check_ksize <- length(unique(signatures_df$ksize))
  check_molecule <- length(unique(signatures_df$molecule))
  check_hash_function <- length(unique(signatures_df$hash_function))
  check_seed <- length(unique(signatures_df$seed))
  check_scaled <- length(unique(signatures_df$seed))
  check_num <- length(unique(signatures_df$seed))

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

#' Check for missing name, NAs in name, or empty name in a data frame of signatures
#'
#' @param signatures_df
#'
#' @return signatures_df with compliant names
#'
#' @examples
#' check_and_edit_names_in_signatures_df(signatures_df)
check_and_edit_names_in_signatures_df <- function(signatures_df){
  # if name isn't in the data frame, or if it is blank, calculate name from basename of filename
  if(!"name" %in% colnames(signatures_df)){
    signatures_df <- signatures_df %>%
      dplyr::mutate(name = basename(filename)) # set names as basename of filename
  }

  # if name is in the data frame, replace missing names (NA) with basename of filename
  if(sum(is.na(signatures_df$name)) > 0){
    signatures_df <- signatures_df %>%
      dplyr::mutate(name = ifelse(is.na(name), basename(filename), name))
  }

  # if name is in data frame, replace it with basename of filename if it contains empty character strings ""
  if("" %in% signatures_df$name){
    signatures_df <- signatures_df %>%
      dplyr::mutate(name = ifelse(name == "", basename(filename), name))
  }

  return(signatures_df)
}

#' Convert a data.frame containing hashes from many sourmash signatures into a upset plot compatible data.frame
#'
#' @description
#' `from_signatures_to_upset_df()` converts a data frame like that produced by binding many data.frames read in by `read_signature()` into a named list of vectors.
#' Each vector in the list is named by the name column in the input signatures data frame and contains the minhashes (mins) for that sample.
#' The `name` column is used for sample names.
#' If the column does not exist, sample names are derived using the base name of filename.
#' If and individual name is blank ("") or NA, that individual name is derived using the base name of filename.
#' The input data frame must contain only one value for each of `ksize`, `scaled`/`num`, `hash_function`, `molecule`, and `seed` as signatures calculated with different values for these parameters are not comparable.
#' The returned data frame will have the minhashes as row names and the shortened filenames as column names and can be plotted with `UpSetR::upset()` or `ComplexUpset::upset()`.
#' If plotting with Complex Upset, additional metadata can be added to the data frame (and therefore the plot) by joining on the values of the minhash rownames.
#'
#' @param signatures_df A data frame of multiple sourmash signatures created by combining many signatures read in by `read_signature()`.
#'
#' @return A an upset plot compliant data frame
#' @export
#'
#' @examples
#' from_signatures_to_upset_df(signatures_df)
from_signatures_to_upset_df <- function(signatures_df){
  # stop if not all parameters to build the sketches are uniform, otherwise the intersections in the upset plot won't make sense
  check_uniform_parameters_in_signatures_df(signatures_df)

  # check and edit signatures names
  signatures_df <- check_and_edit_names_in_signatures_df(signatures_df)

  # turn dataframe into list.
  # Each index in the list is named after the signature it represents.
  # Each index contains a vector of hashes that were in that signature
  upset_list <- signatures_df %>%
    dplyr::group_by(name) %>%
    {setNames(dplyr::group_split(.), dplyr::group_keys(.)[[1]])} # absolutely whack line of code that allows us to set the list names by the group_by variable

  # extract just the minhashes
  upset_list <- lapply(upset_list, function(x) {x$mins})

  # convert list into upset-compliant data.frame
  upset_df <- from_list_to_upset_df(upset_list)

  return(upset_df)
}


#' Visualize the intersection of minhashes in sourmash signatures from many samples
#'
#' @description
#' `plot_signatures_upset` uses `ComplexUpset::upset()` to plot the intersection of minhashes in sourmash signatures from many samples.
#' While based on ggplot2, `ComplexUpset::upset()` has a difficult syntax to parameterize.
#' To produce an alternative visualization, it may be easiest to run `sourmashconsumr::plot_signatures_upset` to retrieve the function source code and alter it.
#'
#' @param upset_df An upset-compliant data frame.
#'
#' @return An upset plot produced by ComplexUpset::upset().
#' @export
#'
#' @examples
#' plot_signatures_upset(upset_df)
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
#' @description
#' `from_signatures_to_rarefaction_df()` rarefies sourmash signatures by sample to assess sequencing depth.
#' The output is a tidy data frame that can be used in ggplot2 plots.
#' It uses `vegan::rarecurve()` to calculate a rarefaction curve for each sample in the input sourmash signatures data frame.
#' The `name` column is used for sample names.
#' If the column does not exist or blank, sample names are derived using the base name of filename.
#' If and individual name is blank ("") or NA, that individual name is derived using the base name of filename.
#' The rarefaction curves are evaluated using the interval of step sample sizes, always including 1 and total sample size.
#' This function is intended to be used on signatures built from read data.
#' It uses abundance information as part of the rarefaction curve calculation.
#' In case signatures were built from reads that have not been k-mer trimmed, there is a filtering step that removes minhashes that are only observed in one sample at abundance 1 as these are likely sequencing errors.
#' This filtering process may invalidate downstream rarefaction curve convergence estimation, as many of these methods evaluate singletons and doubletons in the data set.
#'
#' @param signatures_df A data frame of multiple sourmash signatures created by combining many signatures read in by read_signature(). The data frame must contain the `abundances` column (generated by using the `abund` parameter with `sourmash sketch`).
#' @param step Integer. The step size for samples in rarefaction curve calculation.
#'
#' @return A tidy data frame.
#' @export
#'
#' @examples
#' from_signatures_to_rarefaction_df(signatures_df)
from_signatures_to_rarefaction_df <- function(signatures_df, step = 1){
  # check that the signatures data frame has a column for abundances
  if("abundances" %in% colnames(signatures_df) == FALSE){
    stop("The signatures data.frame does not have an 'abundances' columns.
         Abundances are required to calculate rarefaction cureves.
         If you have abundances in your signatures data.frame, make sure the column is named 'abundances'.")
  }

  # stop if not all parameters to build the sketches are uniform, otherwise the rarefaction curves won't make sense
  check_uniform_parameters_in_signatures_df(signatures_df)

  # format the signatures into a rarecurve-compliant data frame
  signatures_df <- signatures_df %>%
    dplyr::select(name, mins, abundances) %>% # select columns with relevant data
    tidyr::pivot_wider(id_cols = name, values_from = "abundances", names_from = "mins") %>% # transform to wide format
    tibble::column_to_rownames("name") %>%
    replace(is.na(.), 0) # back fill NAs with 0s

  # remove hashes (columns) that are only observed once and in one sample
  # since these rarefaction curves should only be built on reads, these hashes are probably errors
  signatures_df <- signatures_df %>%
    dplyr::select(which(!colSums(signatures_df, na.rm=TRUE) %in% 1))

  # return as a tidy data frame with per-sample information
  # this data frame will include the information needed to plot the rarefaction curves
  rarecurve_df <- vegan::rarecurve(signatures_df, step = step, tidy = T)

  # make the column names more meaningful
  colnames(rarecurve_df) <- c("name", "num_kmers_sampled", "num_kmers_observed")

  return(rarecurve_df)
}

#' Plot a rarefaction curve to determine sequence saturation
#'
#' @description
#' `plot_signatures_rarefaction()` produces a ggplot2 plot from the tidy data frame output by `from_signatures_to_rarefaction_df()`.
#' Since the rarefaction curves may contain many dense points, the `fraction_of_points_to_plot` argument reduces the total number of points rendered which decreases the functions run time.
#' Rarefaction curves that end with a shallow slow likely represent samples that were sequenced deeply enough to capture the k-mer diversity in the sequence.
#'
#' @param rarefaction_df A tidy data frame returned by `from_signatures_to_rarefaction_df()`. Must contain the columns `num_kmers_sampled`, `num_kmers_observed`.
#' @param fraction_of_points_to_plot Integer. 1 every `fraction_of_points_to_plot` will be plotted.
#'
#' @return A ggplot2 object.
#' @export
#'
#' @examples
#' plot_signatures_rarefaction(rarefaction_df)
plot_signatures_rarefaction <- function(rarefaction_df, fraction_of_points_to_plot = 500){
  # check that fraction_of_points_to_plot is >=1
  if(!fraction_of_points_to_plot >= 1){
    stop("Only integers >=1 excepted. If you don't want any filtering, put 1.")
  }

  # we don't need to plot all of the points on the line to get an accurate picture of sequencing saturation
  # This filters to 1/fraction_of_points_to_plot at the correct intervals to maintain curve shape
  # This filtering code will break if given the value 1
  if(fraction_of_points_to_plot != 1){
    rarefaction_df <- rarefaction_df %>%
      dplyr::filter(num_kmers_sampled %% fraction_of_points_to_plot == 1)
  }

  plt <- ggplot2::ggplot(rarefaction_df,
                ggplot2::aes(x = num_kmers_sampled, y = num_kmers_observed)) +
    ggplot2::geom_point() +
    ggplot2::labs(x = "number of sampled k-mers", y = "number of distinct k-mers observed")

  return(plt)
}

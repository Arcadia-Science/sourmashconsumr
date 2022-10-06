#' read_compare_csv
#'
#' @param file path to file output by sourmash compare using the with the --csv flag
#' @param sample_to_rownames Boolean indicating whether sample names should be added to the tibble as a column or a rowname
#' @param ... Arguments passed to readr::read_csv()
#'
#' @return A tibble.
#' @export
#'
#' @examples
#' read_compare_csv("")
read_compare_csv <- function(file, sample_to_rownames = F, ...){
  # guard the use of required packages
  if (!requireNamespace("readr", quietly = TRUE)) {
    stop(
      "Package \"readr\" must be installed to use this function.",
      call. = FALSE
    )
  }

  if (!requireNamespace("dplyr", quietly = TRUE)) {
    stop(
      "Package \"dplyr\" must be installed to use this function.",
      call. = FALSE
    )
  }

  if(sample_to_rownames == F){
  compare_df <- readr::read_csv(file, ...) %>%
    dplyr::mutate(sample = colnames(.), .before = everything())
  } else if(sample_to_rownames == T){

    if (!requireNamespace("tibble", quietly = TRUE)) {
      stop(
        "Package \"tibble\" must be installed to make sample names into row names. Please either install \"tibble\" or set sample_to_rownames to FALSE.",
        call. = FALSE
      )
    }

    compare_df <- readr::read_csv(file, ...) %>%
      dplyr::mutate(sample = colnames(.), .before = everything()) %>%
      tibble::column_to_rownames("sample")
  }
  return(compare_df)
}

#' read gather
#'
#' @param file Path to CSV file or files output by sourmash gather.
#' @param intersect_bp_threshold Integer. Gather matches must have an intersect_bp greater than or equal to this value.
#' @param ... Arguments passed to read_csv().
#'
#' @return A tibble.
#' @export
#'
#' @examples
#' read_gather()
read_gather <- function(file, intersect_bp_threshold, ...){
  if (!requireNamespace("readr", quietly = TRUE)) {
    stop(
      "Package \"readr\" must be installed to use this function.",
      call. = FALSE
    )
  }

  if (!requireNamespace("dplyr", quietly = TRUE)) {
    stop(
      "Package \"dplyr\" must be installed to use this function.",
      call. = FALSE
    )
  }

  if(length(files) > 1){

    if (!requireNamespace("purrr", quietly = TRUE)) {
      stop(
        "Package \"purrr\" must be installed to read multiple files at a time. Please either install \"purrr\" or read one file at a time.",
        call. = FALSE
      )
    }

    # allow the function to read multiple files at once
    taxonomy_annotate_df <- file %>%
      purrr::map_dfr(readr::read_csv, col_types = "ddddddddcccddddcccddcddlddddl", ...) %>%
      dplyr::filter(intersect_bp >= intersect_bp_threshold) %>%
      dplyr::mutate(genome_accession = gsub(" .*", "", name) , .after = "name")

  } else if(length(files) == 1){
    taxonomy_annotate_df <- readr::read_csv(file, col_types = "ddddddddcccddddcccddcddlddddl", ...) %>%
      dplyr::filter(intersect_bp >= intersect_bp_threshold) %>%
      dplyr::mutate(genome_accession = gsub(" .*", "", name) , .after = "name")
  }
}


#' read_taxonomy_annotate
#'
#' @param file Path to CSV file or files output by sourmash taxonomy annotate.
#' @param intersect_bp_threshold Integer. Gather matches must have an intersect_bp greater than or equal to this value.
#' @param separate_lineage Boolean. Read in lineage as a single column or separate each taxonomic level to its own column.
#' @param ... Arguments passed to read_csv().
#'
#' @return A tibble.
#' @export
#'
#' @examples
#' read_taxonomy_annotate()
read_taxonomy_annotate <- function(file, intersect_bp_threshold = 50000, separate_lineage = T, ...){
  if (!requireNamespace("readr", quietly = TRUE)) {
    stop(
      "Package \"readr\" must be installed to use this function.",
      call. = FALSE
    )
  }

  if (!requireNamespace("dplyr", quietly = TRUE)) {
    stop(
      "Package \"dplyr\" must be installed to use this function.",
      call. = FALSE
    )
  }

  if(length(files) > 1){

    if (!requireNamespace("purrr", quietly = TRUE)) {
      stop(
        "Package \"purrr\" must be installed to read multiple files at a time. Please either install \"purrr\" or read one file at a time.",
        call. = FALSE
      )
    }

    # allow the function to read multiple files at once
    taxonomy_annotate_df <- file %>%
      purrr::map_dfr(readr::read_csv, col_types = "ddddddddcccddddcccddcddlddddlc", ...) %>%
      dplyr::filter(intersect_bp >= intersect_bp_threshold) %>%
      dplyr::mutate(n_unique_kmers = (unique_intersect_bp / scaled) * average_abund) %>% # calculate the number of uniquely matched k-mers
      dplyr::mutate(genome_accession = gsub(" .*", "", name) , .after = "name")
  } else if(length(files) == 1){
    taxonomy_annotate_df <- readr::read_csv(file, col_types = "ddddddddcccddddcccddcddlddddlc", ...) %>%
      dplyr::filter(intersect_bp >= intersect_bp_threshold) %>%
      dplyr::mutate(n_unique_kmers = (unique_intersect_bp / scaled) * average_abund) %>% # calculate the number of uniquely matched k-mers
      dplyr::mutate(genome_accession = gsub(" .*", "", name) , .after = "name")
  }
  if(separate_lineage == T){

    if (!requireNamespace("tidyr", quietly = TRUE)) {
      stop(
        "Package \"tidyr\" must be installed to separate taxonomic lineages. Please either install \"tibble\" or set separate_lineage to FALSE.",
        call. = FALSE
      )
    }

    taxonomy_annotate_df <- taxonomy_annotate_df %>%
      tidyr::separate(lineage, into = c("domain", "phylum", "class", "order", "family", "genus", "species", "strain"), sep = ";", remove = F, fill = "right")
  }
  return(taxonomy_annotate_df)
}


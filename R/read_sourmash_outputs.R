#' Read a CSV file output by sourmash compare
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
  if(sample_to_rownames == F){
    compare_df <- readr::read_csv(file, ...) %>%
      dplyr::mutate(sample = colnames(.), .before = everything())
  } else if(sample_to_rownames == T){
    compare_df <- readr::read_csv(file, ...) %>%
      dplyr::mutate(sample = colnames(.), .before = everything()) %>%
      tibble::column_to_rownames("sample")
  }
  return(compare_df)
}

#' Read CSV file or files output by sourmash gather
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
  if(length(file) > 1){
    # allow the function to read multiple files at once
    taxonomy_annotate_df <- file %>%
      purrr::map_dfr(readr::read_csv, col_types = "ddddddddcccddddcccddcddlddddl", ...) %>%
      dplyr::filter(intersect_bp >= intersect_bp_threshold) %>%
      dplyr::mutate(genome_accession = gsub(" .*", "", name) , .after = "name")

  } else if(length(file) == 1){
    taxonomy_annotate_df <- readr::read_csv(file, col_types = "ddddddddcccddddcccddcddlddddl", ...) %>%
      dplyr::filter(intersect_bp >= intersect_bp_threshold) %>%
      dplyr::mutate(genome_accession = gsub(" .*", "", name) , .after = "name")
  }
}


#' Read CSV file or files output by sourmash taxonomy annotate
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
  if(length(file) > 1){
    # allow the function to read multiple files at once
    taxonomy_annotate_df <- file %>%
      purrr::map_dfr(readr::read_csv, col_types = "ddddddddcccddddcccddcddlddddlc", ...) %>%
      dplyr::filter(intersect_bp >= intersect_bp_threshold) %>%
      dplyr::mutate(n_unique_kmers = (unique_intersect_bp / scaled) * average_abund) %>% # calculate the number of uniquely matched k-mers
      dplyr::mutate(genome_accession = gsub(" .*", "", name) , .after = "name")
  } else if(length(file) == 1){
    taxonomy_annotate_df <- readr::read_csv(file, col_types = "ddddddddcccddddcccddcddlddddlc", ...) %>%
      dplyr::filter(intersect_bp >= intersect_bp_threshold) %>%
      dplyr::mutate(n_unique_kmers = (unique_intersect_bp / scaled) * average_abund) %>% # calculate the number of uniquely matched k-mers
      dplyr::mutate(genome_accession = gsub(" .*", "", name) , .after = "name")
  }
  if(separate_lineage == T){
    taxonomy_annotate_df <- taxonomy_annotate_df %>%
      tidyr::separate(lineage, into = c("domain", "phylum", "class", "order", "family", "genus", "species", "strain"), sep = ";", remove = F, fill = "right")
  }
  return(taxonomy_annotate_df)
}

#' Calculate the scaled value for a signature from the max_hash value
#'
#' @param sig_max_hash Max hash value in a signature.
#'
#' @return Integer; a scaled value
#'
#' @examples
get_scaled_for_max_hash <- function(sig_max_hash){
  # sourmash uses the 64-bit hash space of MurmurHash only
  # this is 2 ** 64 - 1 in hexadecimal
  # https://github.com/sourmash-bio/sourmash/blob/a8bd648f60d4e73f3f6bc55cbe70bd174da4a399/src/sourmash/minhash.py#L40
  MINHASH_MAX_HASH = 0xFFFFFFFFFFFFFFFF
  scaled <- round( MINHASH_MAX_HASH / sig_max_hash, digits = 0)
  scaled
}

#' Read sourmash signatures into a dataframe
#'
#' @param file Path to signature (json) file or files output by sourmash sketch (previously sourmash compute).
#' @param compliant Boolean indicating whether signature columns should be compliant; the json fields changed across versions of sourmash. This may drop deprecated columns like 'type' but will allow you to bind many signatures into a single data frame even if they were sketched with different versions of sourmash.
#'
#' @return A tibble.
#' @export
#'
#' @examples
read_signature <- function(file, compliant = TRUE){
  stopifnot(file.exists(file)) # stop if the file doesn't exist
  # read in the signature json file to a dataframe and calculcate the scaled value
  sig_df <- jsonlite::fromJSON(file) %>%
    tidyr::unnest(tidyselect::one_of("signatures")) %>%
    tidyr::unnest(tidyselect::any_of(c("mins", "abundances"))) # any_of allows abundances to be present in signature or not

  # if max_hash is 0, then num was set
  # if num was set, scaled does not need to be calculated
  # max_hash should be the same between all sketches in a sig; I don't think they can be calculated any other way
  if(unique(sig_df$max_hash) != 0){
    sig_df <- sig_df %>%
      dplyr::mutate(scaled = get_scaled_for_max_hash(max_hash)) # calculate the scaled value from max_hash
  }

  if(compliant == TRUE) {
    # define columns that should be selected if compliant = TRUE
    # these are accurate as of sourmash version 4.
    sourmashv4_sig_fields <- c('class', 'email', 'hash_function', 'filename', 'name',
                               'license', 'num', 'ksize', 'seed', 'max_hash', 'mins',
                               'md5sum', 'abundances', 'molecule', 'version')
    sig_df <- sig_df %>%
      dplyr::select(tidyselect::any_of(sourmashv4_sig_fields))
  }

  return(sig_df)
}

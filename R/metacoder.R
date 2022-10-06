

#read_sourmash_taxonomy_annotate()
#pivot_sourmash_taxonomy_annotate_wider()
#parse_tax_data()
#calc_taxon_abund()
#calc_n_samples()



# change sourmash_taxonomy_results to wide format
#' Title
#'
#' @param taxonomy_annotate_tibble
#'
#' @return A tibble in wide format. Variables expanded is n_unique_kmers and colnames are query_name.
#'
#' @examples
#' pivot_wider_taxonomy_annoate(taxonomy_annotate_tibble)
pivot_wider_taxonomy_annotate <- function(taxonomy_annotate_tibble){
  taxonomy_annotate_tibble_wide <- taxonomy_annotate_tibble %>%
    dplyr::select(name, lineage, query_name, n_unique_kmers) %>%
    tidyr::pivot_wider(id_cols = name:lineage, names_from = query_name, values_from = n_unique_kmers)
  return(taxonomy_annotate_tibble_wide)
}

#' taxonomy_annotate_to_metacoder
#'
#' @param taxonomy_annotate_tibble Tibble containing outputs from sourmash taxonomy annotate. If specified, file is ignored. Can contain results from one or many runs of sourmash taxonomy annotate.
#' @inheritParams read_taxonomy_annotate
#' @param summary_level Character.
#' @param ... Arguments passed to metacoder::parse_tax_data().
#' @param groups A tibble or data.frame with distinct query_name values from taxonomy_annotate_tibble in the first column and query groups in the second column.
#' @param groups_prefix Character. Ignored if groups is defined. Used to prefix query_name values for the metacoder function calc_n_samples().
#'
#' @return A metacoder taxmap object.
#' @export
#'
#' @examples
#' taxonomy_annotate_to_metacoder()
taxonomy_annotate_to_metacoder <- function(taxonomy_annotate_tibble = NULL,
                                           file = NULL,
                                           intersect_bp_threshold = 50000,
                                           summary_level = c(NULL, "domain", "phylum", "class", "order", "family", "genus", "species"),
                                           groups = NULL,
                                           groups_prefix = "x",
                                           ...){
  # guard against missing packages
  if (!requireNamespace("dplyr", quietly = TRUE)) {
    stop(
      "Package \"dplyr\" must be installed to use this function.",
      call. = FALSE
    )
  }

  if (!requireNamespace("metacoder", quietly = TRUE)) {
    stop(
      "Package \"metacoder\" must be installed to use this function.",
      call. = FALSE
    )
  }

  # either take in tibble from read_taxonomy_annotate or read in sourmash taxonomy annotate output file(s) directly
  if(missing(taxonomy_annotate_tibble) & missing(file)){
    stop("Neither taxonomy_annotate_tibble or file were specified. Please specify either taxonomy_annotate_tibble or file and retry.")
  }

  # if file is defined but taxonomy_annotate_tibble is missing, read in the files using read_taxonomy_annotate()
  if(missing(taxonomy_annotate_tibble) & !missing(file)){
    taxonomy_annotate_tibble <- read_taxonomy_annotate(file = file, intersect_bp_threshold = intersect_bp_threshold, separate_lineage = F)
  }

  # transform taxonomy annotate tibble into wide format
  taxonomy_annotate_tibble_wide <- pivot_wider_taxonomy_annotate(taxonomy_annotate_tibble)

  # create the metacoder object
  # use the database column to determine which database was used during sourmash taxonomy annotate
  database <- ifelse(grepl(pattern = "genbank", taxonomy_annotate_tibble$filename[1]), "genbank",
                     ifelse(grepl(pattern = "gtdb", taxonomy_annotate_tibble$filename[1]), "gtdb", "userdb"))

  # need this to be able to be user-specified for class_regex and parse_tax_data() in case someone used their own database
  # database <- ()
  if(database == "genbank"){
    metacoder_obj <- metacoder::parse_tax_data(taxonomy_annotate_tibble_wide,
                              class_cols = "lineage", # the column that contains taxonomic information
                              class_sep = ";", # The character used to separate taxa in the classification
                              class_regex = "(.*)", # Regex identifying where the data for each taxon is
                              class_key = c(tax_name = "taxon_name")) # A key describing each regex capture group
  } else if(database == "gtdb"){
    metacoder_obj <- metacoder::parse_tax_data(taxonomy_annotate_tibble_wide,
                              class_cols = "lineage",  # the column that contains taxonomic information
                              class_sep = ";", # The character used to separate taxa in the classification
                              class_regex = "^(.+)__(.+)$", # Regex identifying where the data for each taxon is
                              class_key = c(tax_rank = "info", # A key describing each regex capture group
                                            tax_name = "taxon_name"))
  } else if(!missing(class_regex) & !missing(class_key)){
    metacoder_obj <- metacoder::parse_tax_data(taxonomy_annotate_tibble_wide,
                                               class_cols = "lineage",  # the column that contains taxonomic information
                                               class_sep = ";", # The character used to separate taxa in the classification
                                               class_regex = class_regex,
                                               class_key = class_key)
  } else {
    stop("Could not determine database and no class_regex and class_key defined.\n
         If you used a pre-built sourmash GenBank or GTDB database, please don't alter the values in the database column before running this function.\n
         If you created your own database, please specify class_regex and class_key.")
  }

  # Add information to the metacoder object using metacoder functions

  # set number of unique k-mers assigned to a genome as the abundance information
  metacoder_obj$data$tax_abund <- metacoder::calc_taxon_abund(metacoder_obj, "tax_data", cols = unique(taxonomy_annotate_tibble$query_name))

  # calc_n_samples() calculates the number of samples that contained a taxonomic lineage,
  # and propagates that information up the lineage.
  # we can supply it with sample and group information.
  # the group information is used for adding color to the plots below.

  # by default, samples won't have groups because this information is not encoded in the sourmash taxonomy annotate output.
  # If no groups data.frame is defined by the user, this function uses the query_name as our default groups value.
  # However, we need to modify the query_name value so that metacoder doesn't get confused between a group and a sample.
  if(is.null(groups)){
    groups <- paste0(groups_prefix, colnames(metacoder_obj$data$tax_data)[4:ncol(metacoder_obj$data$tax_data)])
  } else {
    # if user defined a groups dataframe, make sure it's in the same order as samples in the metacoder object
    groups <- groups[order(match(groups[[1]], colnames(metacoder_obj$data$tax_data)[4:ncol(metacoder_obj$data$tax_data)])), ]
    # double check that the groups match with the samples as ordered in the metacoder object
    stopifnot(all.equal(groups[[1]], colnames(metacoder_obj$data$tax_data)[4:ncol(metacoder_obj$data$tax_data)]))
  }

  metacoder_obj$data$tax_occ <- metacoder::calc_n_samples(metacoder_obj, "tax_abund",
                                                          groups =  groups[[2]],
                                                          cols = unique(taxonomy_annotate_tibble$query_name))

  return(metacoder_obj)


}

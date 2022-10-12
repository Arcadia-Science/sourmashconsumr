#' Pivot sourmash_taxonomy_results to wide format
#'
#' @param taxonomy_annotate_tibble
#'
#' @return A tibble in wide format. Variables expanded is n_unique_kmers and colnames are query_name.
#'
#' @examples
#' pivot_wider_taxonomy_annotate(taxonomy_annotate_tibble)
pivot_wider_taxonomy_annotate <- function(taxonomy_annotate_tibble){
  taxonomy_annotate_tibble_wide <- taxonomy_annotate_tibble %>%
    dplyr::select_if(colnames(.) %in% c("genome_accession", "lineage", "query_name", "n_unique_kmers")) %>% # use select_if to allow genome_accession to be missing, as it won't be present in agglomerated columns
    tidyr::pivot_wider(names_from = query_name, values_from = n_unique_kmers) # leverage default behavior to have everything be an id col other than names_from and values_from
  return(taxonomy_annotate_tibble_wide)
}

#' Import the output of sourmash taxonomy annotate into a taxmap metacoder object
#'
#' @param taxonomy_annotate_tibble Tibble containing outputs from sourmash taxonomy annotate. If specified, file is ignored. Can contain results from one or many runs of sourmash taxonomy annotate.
#' @inheritParams read_taxonomy_annotate
#' @param summary_level Character. NULL by default, meaning no summarization is done. Valid options are "domain", "phylum", "class", "order", "family", "genus", and "species". When a valid option is supplied, k-mer counts are agglomerated to that level before metacoder object is created.
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
                                           summary_level = NULL,
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

  # figure out what database was used before manipulating the taxonomy_annotate_tibble:
  # use the database column to determine which database was used during sourmash taxonomy annotate
  database <- ifelse(grepl(pattern = "genbank", taxonomy_annotate_tibble$filename[1]), "genbank",
                     ifelse(grepl(pattern = "gtdb", taxonomy_annotate_tibble$filename[1]), "gtdb", "userdb"))

  # if summary_level is defined, parse lineage and agglomerate counts to that level of taxonomy
  if(!is.null(summary_level)){
    # make sure only except arguments are returned
    if(!summary_level %in% c("domain", "phylum", "class", "order", "family", "genus", "species")){
      stop("Unrecognized string passed to summary_level. Please use one of species, genus, family, order, class, phylum, or domain.")
    }
    # agglomerate to level of taxonomy
    if(summary_level == "domain"){
      agglom_cols <- c("query_name", "domain")
    } else if(summary_level == "phylum"){
      agglom_cols <- c("query_name", "domain", "phylum")
    } else if(summary_level == "class"){
      agglom_cols <- c("query_name", "domain", "phylum", "class")
    } else if(summary_level == "order"){
      agglom_cols <- c("query_name", "domain", "phylum", "class", "order")
    } else if(summary_level == "family"){
      agglom_cols <- c("query_name", "domain", "phylum", "class", "order", "family")
    } else if(summary_level == "genus"){
      agglom_cols <- c("query_name", "domain", "phylum", "class", "order", "family", "genus")
    } else if(summary_level == "species"){
      agglom_cols <- c("query_name", "domain", "phylum", "class", "order", "family", "genus", "species")
    }

    taxonomy_annotate_tibble <- taxonomy_annotate_tibble %>%
      dplyr::select(genome_accession, lineage, query_name, n_unique_kmers) %>%
      tidyr::separate(lineage, into = c("domain", "phylum", "class", "order", "family", "genus", "species", "strain"), sep = ";", remove = F, fill = "right") %>%
      dplyr::group_by_at(dplyr::vars(dplyr::all_of(agglom_cols))) %>%
      dplyr::summarize(n_unique_kmers = sum(n_unique_kmers)) %>%
      dplyr::ungroup() %>%
      tidyr::unite(lineage, all_of(agglom_cols[-1]), sep = ";", remove = TRUE) %>%
      dplyr::select(lineage, query_name, n_unique_kmers)

  }

  # transform taxonomy annotate tibble into wide format
  taxonomy_annotate_tibble_wide <- pivot_wider_taxonomy_annotate(taxonomy_annotate_tibble)

  # create the metacoder object
  # need this to be able to be user-specified for class_regex and parse_tax_data() in case someone used their own database
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

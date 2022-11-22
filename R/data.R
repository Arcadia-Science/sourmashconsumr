#' Gut microbiome sourmash compare data frame
#'
#' A similarity matrix between six stool microbiome shotgun metagenome samples.
#'
#' @format ## `gut_compare_df`
#' A data frame with 6 rows and 7 columns:
#' \describe{
#'   \item{sample}{Sample name}
#'   \item{SRR5935765, SRR5936131, SRR5936197, SRR5946920, SRR5946923, SRR5947006}{Sequence Read Archive run accessions. Aligned values specify similarity between two samples.}
#'   ...
#' }
#' @source <https://github.com/Arcadia-Science/sourmashconsumr/blob/main/data-raw/00_sourmash_commands.sh>
"gut_compare_df"


#' Gut microbiome sourmash gather data frame
#'
#' The results from running sourmash gather on each of six stool microbiome
#' metagenomes against the GTDB rs207 representatives database.
#'
#' @format ## `gut_gather_df`
#' A data frame with 1,062 rows and 30 columns:
#' \describe{
#'   \item{intersect_bp}{Numeric. Estimated number of intersected base pairs between a metagenome and a genome in a database.}
#'   \item{f_orig_query}{Numeric. Fraction of the original query that belongs to the match.}
#'   \item{f_match}{Numeric. Fraction of the matched genome in the leftover query.}
#'   \item{f_unique_to_query}{Numeric. Fraction of the query that uniquely belongs to the match.}
#'   \item{f_unique_weighted}{Numeric. Abundance-weighted fraction of the query that uniquely belongs to the match.}
#'   \item{average_abund}{Numeric. Average abundance of k-mers in the metagenome that were in the match.}
#'   \item{median_abund}{Numeric. Median abundance of k-mers in the metagenome that were in the match.}
#'   \item{std_abund}{Numeric. Standard deviation of abundance of k-mers in the metagenome that were in the match.}
#'   \item{filename}{Character. File path for the database on the computer that sourmash gather was executed on.}
#'   \item{name}{Character. Name of matched genome in the sourmash gather database.}
#'   \item{genome_accession}{Character. Genome accession solved by cutting of the name variable at the first space.}
#'   \item{md5}{Character. MD5 hash for the matched genome sketch.}
#'   \item{f_match_orig}{Numeric. Fraction of the matched genome in the original query prior to gather subtraction.}
#'   \item{unique_intersect_bp}{Numeric. Estimated number of uniquely intersected base pairs between a metagenome and a genome in a database.}
#'   \item{gather_result_rank}{Numeric. Rank of match in gather results.}
#'   \item{remaining_bp}{Numeric. Remaining base pairs in the query after the match is removed.}
#'   \item{query_filename}{Character. File name for the query derived from the query sketch.}
#'   \item{query_name}{Character. Name of the query.}
#'   \item{query_md5}{Character. MD5 hash for the query.}
#'   \item{query_bp}{Character. Number of base pairs in the query.}
#'   \item{ksize}{Character. K-mer size used for sourmash gather.}
#'   \item{moltype}{Character. Molecule type used for sourmash gather.}
#'   \item{scaled}{Numeric. Scaled value that sourmash gather was performed at.}
#'   \item{query_n_hashes}{Numeric. Number of hashes (k-mers) in the query.}
#'   \item{query_abundance}{Logical. Whether hash (k-mer) abundance information was a part of the query sketch.}
#'   \item{query_containment_ani}{Numeric. Containment between the query and the match.}
#'   \item{match_containment_ani}{Numeric. Containment between the match and the query.}
#'   \item{average_containment_ani}{Numeric. Average of the two containment metrics query_containment_ani and match_containment_ani.}
#'   \item{max_containment_ani}{Numeric. Maximum containment between the two containment metrics query_containment_ani and match_containment_ani.}
#'   \item{potential_false_negative}{Logical. Whether the containment estimate is a potential false negative.}
#'   ...
#' }
#' @source <https://github.com/Arcadia-Science/sourmashconsumr/blob/main/data-raw/00_sourmash_commands.sh>
"gut_gather_df"


#' Gut microbiome sourmash taxonomy annotate data frame
#'
#' The results from running sourmash gather on each of six stool microbiome
#' metagenomes against the GTDB rs207 representatives database and then
#' assigning taxonomy using sourmash taxonomy annotate.
#'
#' @format ## `gut_taxonomy_annotate_df`
#' A data frame with 1,062 rows and 40 columns:
#' \describe{
#'   \item{intersect_bp}{Numeric. Estimated number of intersected base pairs between a metagenome and a genome in a database.}
#'   \item{f_orig_query}{Numeric. Fraction of the original query that belongs to the match.}
#'   \item{f_match}{Numeric. Fraction of the matched genome in the leftover query.}
#'   \item{f_unique_to_query}{Numeric. Fraction of the query that uniquely belongs to the match.}
#'   \item{f_unique_weighted}{Numeric. Abundance-weighted fraction of the query that uniquely belongs to the match.}
#'   \item{average_abund}{Numeric. Average abundance of k-mers in the metagenome that were in the match.}
#'   \item{median_abund}{Numeric. Median abundance of k-mers in the metagenome that were in the match.}
#'   \item{std_abund}{Numeric. Standard deviation of abundance of k-mers in the metagenome that were in the match.}
#'   \item{filename}{Character. File path for the database on the computer that sourmash gather was executed on.}
#'   \item{name}{Character. Name of matched genome in the sourmash gather database.}
#'   \item{genome_accession}{Character. Genome accession solved by cutting of the name variable at the first space.}
#'   \item{md5}{Character. MD5 hash for the matched genome sketch.}
#'   \item{f_match_orig}{Numeric. Fraction of the matched genome in the original query prior to gather subtraction.}
#'   \item{unique_intersect_bp}{Numeric. Estimated number of uniquely intersected base pairs between a metagenome and a genome in a database.}
#'   \item{gather_result_rank}{Numeric. Rank of match in gather results.}
#'   \item{remaining_bp}{Numeric. Remaining base pairs in the query after the match is removed.}
#'   \item{query_filename}{Character. File name for the query derived from the query sketch.}
#'   \item{query_name}{Character. Name of the query.}
#'   \item{query_md5}{Character. MD5 hash for the query.}
#'   \item{query_bp}{Character. Number of base pairs in the query.}
#'   \item{ksize}{Character. K-mer size used for sourmash gather.}
#'   \item{moltype}{Character. Molecule type used for sourmash gather.}
#'   \item{scaled}{Numeric. Scaled value that sourmash gather was performed at.}
#'   \item{query_n_hashes}{Numeric. Number of hashes (k-mers) in the query.}
#'   \item{query_abundance}{Logical. Whether hash (k-mer) abundance information was a part of the query sketch.}
#'   \item{query_containment_ani}{Numeric. Containment between the query and the match.}
#'   \item{match_containment_ani}{Numeric. Containment between the match and the query.}
#'   \item{average_containment_ani}{Numeric. Average of the two containment metrics query_containment_ani and match_containment_ani.}
#'   \item{max_containment_ani}{Numeric. Maximum containment between the two containment metrics query_containment_ani and match_containment_ani.}
#'   \item{potential_false_negative}{Logical. Whether the containment estimate is a potential false negative.}
#'   \item{lineage}{Character. Full taxonomic lineage of the gather match. Each taxonomic level is separated by a semicolon.}
#'   \item{domain}{Character. Domain in the taxonomic lineage.}
#'   \item{phylum}{Character. Phylum in the taxonomic lineage.}
#'   \item{class}{Character. Class in the taxonomic lineage.}
#'   \item{order}{Character. Order in the taxonomic lineage.}
#'   \item{family}{Character. Family in the taxonomic lineage.}
#'   \item{genus}{Character. Genus in the taxonomic lineage.}
#'   \item{species}{Character. Species in the taxonomic lineage.}
#'   \item{strain}{Character. Strain in the taxonomic lineage.}
#'   \item{n_unique_kmers}{Numeric. Abundance-weighted number of unique k-mers attributable to the gather match.}
#'   ...
#' }
#' @source <https://github.com/Arcadia-Science/sourmashconsumr/blob/main/data-raw/00_sourmash_commands.sh>
"gut_taxonomy_annotate_df"


#' Gut microbiome sourmash signatures data frame
#'
#' Sourmash sketches from six stool microbiome shotgun metagenome samples parsed
#' into a data frame.
#'
#' @format ## `gut_signatures_df`
#' A data frame with 15,339 rows and 16 columns:
#' \describe{
#'   \item{class}{Character. Class of the sketch.}
#'   \item{email}{Character. Email recorded in the sketch.}
#'   \item{hash_function}{Character. Hash function used for the sketch.}
#'   \item{filename}{Character. Name of the file that was sketched.}
#'   \item{name}{Character. Name of the sketch.}
#'   \item{license}{Character. License for the sketch.}
#'   \item{num}{Integer. Number of hashes specified to be in the sketch.}
#'   \item{ksize}{Integer. K-mer size represented in the sketch.}
#'   \item{seed}{Integer. Seed for the hash function.}
#'   \item{max_hash}{Numeric. Maximum hash.}
#'   \item{scaled}{Numeric. Scaled value for the sketch.}
#'   \item{mins}{Numeric. MinHash, or hash, that represents each k-mer sketched from the sample.}
#'   \item{md5sum}{Character. MD5 hash for the sketch.}
#'   \item{abundances}{Integer. Abundance of hash in sample.}
#'   \item{molecule}{Character. Molecule type represented by sketch.}
#'   \item{version}{Numeric. Sketch version.}
#'   ...
#' }
#' @source <https://github.com/Arcadia-Science/sourmashconsumr/blob/main/data-raw/00_sourmash_commands.sh>
"gut_signatures_df"

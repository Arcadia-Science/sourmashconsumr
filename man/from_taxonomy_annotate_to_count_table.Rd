% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/from_taxonomy_annotate_to_phyloseq.R
\name{from_taxonomy_annotate_to_count_table}
\alias{from_taxonomy_annotate_to_count_table}
\title{Transform the output of sourmash taxonomy annotate into a phyloseq OTU table (count table)}
\usage{
from_taxonomy_annotate_to_count_table(taxonomy_annotate_df)
}
\arguments{
\item{taxonomy_annotate_df}{Data frame containing outputs from sourmash taxonomy annotate.}
}
\value{
A phyloseq otu_table object
}
\description{
Transform the output of sourmash taxonomy annotate into a phyloseq OTU table (count table)
}
\examples{
\dontrun{
from_taxonomy_annotate_to_count_table(taxonomy_annotate_df)
}
}

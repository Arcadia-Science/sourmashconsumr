% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot_taxonomy_annotate.R
\name{from_taxonomy_annotate_to_multi_strains}
\alias{from_taxonomy_annotate_to_multi_strains}
\title{Detect whether multiple strains of the same species are present in a sample.}
\usage{
from_taxonomy_annotate_to_multi_strains(
  taxonomy_annotate_df,
  plot_threshold = 0.02
)
}
\arguments{
\item{taxonomy_annotate_df}{Data frame containing outputs from sourmash taxonomy annotate.
Can contain results from one or many runs of sourmash taxonomy annotate.}

\item{plot_threshold}{f_match threshold for plotting a genome match.
This threshold is for plotting only.}
}
\value{
A named list.
The first object in the list `candidate_species_with_multiple_strains` summarizes the query sample and species names that may have an f_match >= 1.1.
The second object `plt` is a ggplot object that summarizes the query samples and species that may contain multiple strains.
The third object `plt_data` contains the data that is used to produce the ggplot object.
}
\description{
`from_taxonomy_annotate_to_multi_strains()` uses the output of sourmash taxonomy annotate to detect whether multiple strains of the same species are present in a sample (e.g. a metagenome).
The function uses the `f_match` metric produced by sourmash gather to predict if multiple strains are present.
`f_match` is the fraction of a matched genome that is contained within the query sample.
This function sums over all `f_match` values for all matched genomes of a given species and detects when we see more genomic segments than we would expect to see if only one strain were truly present.
Most frequently, if multiple genomes form one species are identified by sourmash gather, the single genome in a query sample is different from genomes in databases so multiple genomes that each cover a distinct portion of the genome that is truly present are returned as matches.
However, for each species, the fraction of a single genome that matched against the query decreases with each successive match.
This is consistent with the idea of pangenomes -- the single true genome in our query matches many different parts of genomes in the database.
When this happens, the sum fraction matched (`f_match`) by all genomes within a species within a single query should not exceed ~1.
When we run sourmash gather on single genomes, the real value we observe for the sum f_match ranged between:
0.046 (for genomes that only had ~genus-level relatives in the database) and 1.04 (for genomes that had many close matches in the database).
We feel reasonably confident that it is a good starting place for strain-level analysis (see details below for caveats) but additional validated methods should be used to confirm findings.
}
\details{
In this context, we refer to _strain_ as any sub-species level variation.
This function uses genomes as the unit of strain -- each genome is considered a different strain.
Sourmash gather compares a query (e.g. a metagenome) against a database (such as GenBank microbial genomes) and provides the minimum set of genomes that cover all of the k-mers in the query that are in the database.
At least two things could be happening when sourmash gather returns two genomes of the same species (e.g. different strains) as a match to the same metagenome sample:
1. Both strains may be present in the metagenome
2. Only one strain may be truly present in the metagenome, but that strain is not contained within our current reference database.
Instead, pieces of that strain's genome are in other genomes in the database.
The genome in the database that contains the largest overlap with the genome in the metagenome is returned first as the best match.
Then, other genomes in the database are returned that match other portions of the metagenome strain's genome that wasn't contained in the best match.
In reality, some combination of these two things probably happens.
Keep in mind:
1. We can't detect strain variation if there is only one genome for a given species in the database using sourmash gather/taxonomy alone (variant calling tools or tools that look at variation in assembly graphs would be more successful for this use case).
2. The sourmash gather/taxonomy results alone should not be used to conclusively detect strain variation, but it is a good place to start to figure out where to dig in deeper.
}
\examples{
\dontrun{
from_taxonomy_annotate_to_multi_strains()
}
}

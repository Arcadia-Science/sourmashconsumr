% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot_gather.R
\name{plot_gather_upset}
\alias{plot_gather_upset}
\title{Visualize the intersection of genome accessions in sourmash gather results from many samples}
\usage{
plot_gather_upset(
  upset_df,
  color_by_database = FALSE,
  gather_df = NULL,
  palette = NULL,
  ...
)
}
\arguments{
\item{upset_df}{An upset plot compliant data frame like that produced by from_gather_to_upset_df.}

\item{color_by_database}{Boolean indicating whether to fill the upset plot instersection bar plot by the database the results came from. FALSE by default.}

\item{gather_df}{Gather results passed to from_gather_to_upset_df to create the input data frame. NULL unless color_by_database is set to TRUE.}

\item{palette}{Optional argument specifying the palette to use. Ignored unless color_by_database is set to TRUE. Defaults to RColorBrewer Set2 if color_by_database is set to TRUE and palette is not specified.}

\item{...}{Arguments passed to ComplexUpset::upset().}
}
\value{
A ComplexUpset plot.
}
\description{
`plot_gather_upset` uses `ComplexUpset::upset()` to plot the intersection of genome accessions in sourmash gather results from many samples.
The plot can optionally be colored by the database from which the genome accession originated.
While based on ggplot2, `ComplexUpset::upset()` has a difficult syntax to parameterize.
To produce an alternative visualization, it may be easiest to run `sourmashconsumr::plot_gather_upset` to retrieve the function source code and alter it.
}
\examples{
\dontrun{
plot_gather_upset(upset_df)
}
}

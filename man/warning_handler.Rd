% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/read_sourmash_outputs.R
\name{warning_handler}
\alias{warning_handler}
\title{Warning handler to suppress warning message that contains "The following named parsers don't match the column names:"}
\usage{
warning_handler(w)
}
\arguments{
\item{w}{Warning passed to the function by withCallingHandlers()}
}
\value{
The object created by the first argument passed to the function without a warning if the warning contains "The following named parsers don't match the column names:"
}
\description{
Function to suppress warning messages produced by `read_gather()` and `read_taxonomy_annotate()` when results were generated by earlier versions of sourmash that didn't output all of the columns output by the current version of sourmash gather.
}
\examples{
\dontrun{
withCallingHandlers(df <- read_csv(file, col_types = col_types),
                    warning = warning_handler))
}
}


<!-- README.md is generated from README.Rmd. Please edit that file -->

# sourmashconsumr

<!-- badges: start -->

[![R-CMD-check](https://github.com/Arcadia-Science/sourmashconsumr/actions/workflows/check-standard.yaml/badge.svg)](https://github.com/Arcadia-Science/sourmashconsumr/actions/workflows/check-standard.yaml)
[![Codecov test
coverage](https://codecov.io/gh/Arcadia-Science/sourmashconsumr/branch/main/graph/badge.svg)](https://app.codecov.io/gh/Arcadia-Science/sourmashconsumr?branch=main)
<!-- badges: end -->

The goal of sourmashconsumr is to parse, analyze, and visualize the
outputs of the [sourmash python
package](https://sourmash.readthedocs.io/en/latest/). The
sourmashconsumr package is still under active development.

## Installation

You can install the development version of sourmashconsumr from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("Arcadia-Science/sourmashconsumr")
```

Eventually, we hope to release sourmashconsumr on CRAN and to provide a
bioconda package. We’ll update these instructions once we’ve done that.

## Usage

See the vignette for full instructions on how to run the sourmashconsumr
package (coming soon!).

To access the functions in the sourmashconsumr package, you can load it
with:

``` r
library(sourmashconsumr)
## basic example code
```

The sourmashconsumr package contains a variety of functions to work with
the outputs of the sourmash pythong package. `read*` functions read and
parse the outputs of sourmash into data frames. `plot*` functions
visualize the sourmash outputs. `from_*to_*` functions convert sourmash
data frames into the formats of other popular R packages
(e.g. phyloseq). sourmashconsumr also contains a few utility functions
that work with other datatype (at this time, mostly upset data frames).
See a complete list of functions below, sorted by data type that the
function operates on:

- signatures (output by `sourmash sketch` or `sourmash compute`):
  - `read_signature()`
  - upset plots: `from_signatures_to_upset_df()`,
    `plot_signatures_upset()`
  - rarefaction plots for signatures sketched from reads:
    `from_signatures_to_rarefaction_df()`,
    `plot_signatures_rarefaction()`
- sourmash compare csv:
  - `read_compare()`
  - MDS plot: `make_compare_mds()`, `plot_compare_mds()`
  - heatmap: `plot_compare_heatmap()`
- sourmash taxonomy annotate csv:
  - `read_taxonomy_annotate()`
  - taxonomy agglomeration: `tax_glom_taxonomy_annotate()`
  - upset plot: `from_taxonomy_annotate_to_upset_inputs()`,
    `plot_taxonomy_annotate_upset()`
  - sankey plot: `plot_taxonomy_annotate_sankey()`
  - time series alluvial plot: `plot_taxonomy_annotate_ts_alluvial()`
  - to metacoder: `from_taxonomy_annotate_to_metacoder()`
  - to phyloseq: `from_taxonomy_annotate_to_phyloseq()`
- sourmash gather csv:
  - `read_gather()`
  - barchart: `plot_gather_classified()`
  - upset plot: `from_gather_to_upset_df()`, `plot_gather_upset()`
- upset utilities:
  - `from_list_to_upset_df()`
  - `from_upset_df_to_intersection_members()`
  - `from_upset_df_to_intersection_summary()`
  - `from_upset_df_to_intersections()`

## Citation

sourmashconsumr doesn’t have a citation yet, but sourmash does! If you
use sourmash in your work, please cite: [Brown and Irber (2016),
doi:10.21105/joss.00027.](https://joss.theoj.org/papers/10.21105/joss.00027)

If you’d like more information on how sourmash works, please see the
following publications: \* For a general background on how sourmash
works and examples of how to use it: [Large-scale sequence comparisons
with sourmash](https://f1000research.com/articles/8-1006/v1) \* For a
mathematical description of FracMinHash and a demonstration of the
accuracy of sourmash gather: [Lightweight compositional analysis of
metagenomes with FracMinHash and minimum metagenome
covers](https://www.biorxiv.org/content/10.1101/2022.01.11.475838v2.abstract)

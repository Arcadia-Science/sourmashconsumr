
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
# install.packages("remotes")
remotes::install_github("Arcadia-Science/sourmashconsumr")
```

Eventually, we hope to release sourmashconsumr on CRAN and to provide a
conda-forge package. We’ll update these instructions once we’ve done
that.

## Usage

See the vignette for full instructions on how to run the sourmashconsumr
package (coming soon!).

To access the functions in the sourmashconsumr package, you can load it
with:

``` r
library(sourmashconsumr)
```

The sourmashconsumr package contains a variety of functions to work with
the outputs of the sourmash python package. The table below summarizes
which sourmash outputs the sourmashconsumr package operates on and the
functions that are available. For a complete list of functions in the
sourmashconsumr package, see the
[documentation](https://arcadia-science.github.io/sourmashconsumr/reference/index.html).

<img src="https://i.imgur.com/UfuiAhw.png" width="750" />

## Developer documentation

The sourmashconsumr package follows package developer conventions laid
out in <https://r-pkgs.org/>, and changes can be contributed to the code
base using pull requests. For more information on how to contribute, see
the [developer documentation](devdoc.md).

## Citation

- If you use sourmashconsumr in your work, please cite [DOI:
  10.57844/arcadia-1896-ke33](https://arcadia-research.pubpub.org/pub/resource-sourmashconsumr).
- If you use sourmash in your work, please cite [DOI:
  10.21105/joss.00027](https://joss.theoj.org/papers/10.21105/joss.00027).

If you’d like more information on how sourmash works, please see the
following publications:

- For a general background on how sourmash works and examples of how to
  use it: [Large-scale sequence comparisons with
  sourmash](https://f1000research.com/articles/8-1006/v1)
- For a mathematical description of FracMinHash and a demonstration of
  the accuracy of sourmash gather: [Lightweight compositional analysis
  of metagenomes with FracMinHash and minimum metagenome
  covers](https://www.biorxiv.org/content/10.1101/2022.01.11.475838v2.abstract)

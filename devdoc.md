# Developer documentation

The sourmashconsumr package follows package developer conventions laid out in [https://r-pkgs.org/](https://r-pkgs.org/).
It uses this GitHub repository, the software for R and RStudio, and the R packages `devtools`, `testthat`, and `usethis` for development. 
The packages in the `DESCRIPTION` file are dependencies.

## Obtaining the sourmashconsumr source code and installing package and development dependencies

### Obtaining the source code

+ Developers who are not members of the [Arcadia-Science](https://github.com/Arcadia-Science/) organization: fork the repo to your own GitHub user account. Clone the repository locally.
+ Developers who are members of the [Arcadia-Science](https://github.com/Arcadia-Science/) organization: clone the repository and create a branch.
Clone the repository locally using `git clone`.

See [contributing section](#contributing-changes-back-to-the-sourmashconsumr-repository) below for instructions on how to contribute code changes back to the sourmashconsumr repository.

### Installing R and RStudio

Make sure you have R and RStudio installed (see [here](https://rstudio-education.github.io/hopr/starting.html) for operating-system specific instructions).

### Opening the sourmashconsumr R project

You'll then need to open the R project for the `sourmashconsumr` project. 
You can do this by double clicking the `sourmashconsumr.Rproj` file in your file finder or by using `File` -> `Open Project...` and then finding and opening the `sourmashconsumr.Rproj` file.

### Installing development dependencies

To start developing the package, you'll need to install `devtools`, `testthat`, and `usethis`. 

You can use the following R command to install these packages:

```
install.packages(c("devtools", "testthat", "usethis"))
```

To be able to run all functions, and to run all of the tests, you will also need to install the package dependencies. 
See the instructions below for options for how to install these dependencies.

### Installing the package dependencies

The sourmashconsumr package has a lot of dependencies as it tries to make the outputs of the sourmash python package interoperable with biological computing packages already in the R ecosystem. You can install these dependencies however you like; they're all documented in the `DESCRIPTION` file in the root of this repository. Below we describe different strategies for installing the dependencies.

#### Installing package dependencies using R functions

Installing the sourmashconsumr package from GitHub will trigger missing dependencies to be installed.

```
install.packages("devtools")
devtools::install_github("Arcadia-Science/sourmashconsumr")
```

Some of the package dependencies are available on Bioconductor (namely phyloseq), so if installation of dependencies fails through `install_github`, you may have more success using `BiocManager`:
```
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("Arcadia-Science/sourmashconsumr")
```

#### Managing the development environment with conda

If the above methods fail, conda can be used to install the dependencies.
You can find operating system-specific instructions for installing miniconda [here](https://docs.conda.io/en/latest/miniconda.html).
Once you have conda and [mamba](https://mamba.readthedocs.io/en/latest/) installed, you can create an environment with all of the sourmashconsumr dependencies and developer packages installed using the following bash command:

```
mamba env create -n sourmashconsumr --file devenv.yml
```

Then, activate the environment:

```
conda activate sourmashconsumr
```

You will still need to install RStudio after creating the conda environment. 
Then, if you're using a Mac, you can open the sourmashconsumr Rproject by running `open sourmashconsumr.Rproj` from the activated environment you just created.

```
open sourmashconsumr.Rproj
```

On a Linux or Windows machine, running `rstudio` from the activated conda environment may open RStudio (untested).

```
rstudio
```

## Developing 

sourmashconsumr is developed on GitHub at [Arcadia-Science/sourmashconsumr](https://github.com/Arcadia-Science/sourmashconsumr/).
Changes can be contributed to sourmashconsumr via [pull requests](https://github.com/Arcadia-Science/sourmashconsumr/pulls).

### Getting started

To start developing, open the sourmashconsumr.Rproj in RStudio (see above).
Then, load the `devtools` package and use the `devtools::load_all()` function to load sourmashconsumr.

```
library(devtools)
load_all()
```

You can then open the file you want to change and make changes.

### Adding or changing functions

Core package functions can be found in the `R` folder. 
Functions are grouped by action (e.g. `read_sourmash_outputs.R` contains all of the `read_*()` functions) or by data type (e.g. functions that involve plotting the output of sourmash taxonomy annotate are in the `plot_taxonomy_annotate.R` file).

#### Naming functions

The functions in the sourmashconsumr package follow a naming scheme
Functions that are exported (e.g. user-facing) are named by the action completed by the function, the sourmash output type the act on, and if relevant, a description of the action taken.

+ **Action words**:
    + `read`
    + `plot`
    + `from`
+ **sourmash output types**:
    + `signature`
    + `compare_csv`
    + `gather`
    + `taxonomy_annotate`
+ **example actions**:
    + `to_metacoder`
    + `upset`
    + `heatmap`
    + `mds`
    
Functions that are not exported do not follow a naming scheme but strive to be fully descriptive of their actions, and when possible use the **sourmash output types** to make it clear what type of data the internal function operates on. 

+ **examples of internal functions**
   + `check_compare_df_sample_col_and_move_to_rowname()`
   + `check_and_edit_names_in_signatures_df()`
   + `check_uniform_parameters_in_signatures_df()`
   + `make_agglom_cols()`
   + `make_expression()`
   + `get_scaled_for_max_hash()`
   + `pivot_wider_taxonomy_annotate()`

#### Other things

For the `plot*` functions, sourmashconsumr tries not to put anything in the `theme()` layer so that users can control the look of their plot with out conflicts.
This isn't always do-able, but we strive to follow this rule.

Layers that label the plot should be controlled to by a boolean argument `label` so that users can turn off labels and re-add them themselves so they can control the look of those labels (italicized, size, etc.).
When there is a label boolean as a function argument, the `@param` string documenting that parameter should contain an example of how to add the default labels back to the plot as some of these are esoteric and would be hard to divine without an example.

### Tests

The sourmashconsumr package uses unit tests to make sure that code changes don't break existing functions.
Tests and test files are located in the `tests/testthat` folder.
To run all tests, you can use the `testthat::test()` function:

```
test()
```

Alternatively, you can run tests related to the active file:
```
test_active_file()
```

### Documentation

Documentation is coordinated by roxygen and the `devtools::document()` function.
The documentation for each function is contained in roxygen blocks (`#'`) above the function.
These blocks are automatically parsed by the `document()` function to produce the files in the `man/*` directory. 
Therefore, the files in the `man/*` directory should not be manually edited.

The only exception to this is the packages listed in the DESCRIPTION file.
These packages are added to the DESCRIPTION file using the `usethis::use_package()` function.

### Vignette

The vignette html file is built from the Rmd file `vignettes/sourmashconsumr.Rmd`.
The vignette can be built using the following command:

```
devtools::build_rmd("vignettes/sourmashconsumr.Rmd")
```

If you make changes to the vignette, make sure you build it locally and that your changes appear how you want them to before you push changes to GitHub.
The html file should not be pushed to GitHub.

### GitHub README

The README.md is created from a README.Rmd file.
To edit the README.md file, make changes to the README.Rmd file and knit it.

### pkgdown package documentation site

sourmashconsumr uses [pkgdown](https://pkgdown.r-lib.org/index.html) to build the sourmashconsumr documentation.
To build the package documentation site locally, you can use:
```
usethis::use_pkgdown()
pkgdown::build_site()
```

You shouldn't have to do this next part, but to set up GitHub Actions so the site would be automatically built, we used the following code:
```
# set up git credentials
usethis::gh_token_help()
usethis::create_github_token()
gitcreds::gitcreds_set()
# use gh-pages and github actions to build site
usethis::use_pkgdown_github_pages()
```

### Checking all the things

After making changes, you should run `devtools::check()` locally to make sure all tests pass, documentation is updated, and there are no other failures or warning in the local build process.

```
check()
```

## Contributing changes back to the sourmashconsumr repository

Contributions can be made to the sourmashconsumr package using [pull requests](https://github.com/Arcadia-Science/sourmashconsumr/pulls).
While not required, you may want to open an [issue](https://github.com/Arcadia-Science/sourmashconsumr/issues) first to discuss your proposed changes.
Pull requests should be opened against the `main` branch.
When you open a pull request from your branch or fork, GitHub Actions will automatically launch continuous integration checks to make sure that the proposed changes pass R-CMD check across different operating systems (Mac, Windows, Linux).
All checks must pass for new code to be merged into the `main` branch. 
The CI actions will also calculate changes to code coverage, or how much of the code base is covered by tests.
All exported functions should have a corresponding test.

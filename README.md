# mia - Microbiome analysis <img src="man/figures/mia_logo.png" align="right" width="120" />

<!-- badges: start -->

[![R-CMD-check-bioc-devel](https://github.com/microbiome/mia/workflows/R-CMD-check-bioc-devel/badge.svg)](https://github.com/microbiome/mia/actions)
<!--[![R-CMD-check-bioc](https://github.com/microbiome/mia/actions/workflows/check-bioc.yml/badge.svg)](https://github.com/microbiome/mia/actions/workflows/check-bioc.yml)-->
[![Codecov test
coverage](https://codecov.io/gh/microbiome/mia/branch/master/graph/badge.svg)](https://codecov.io/gh/microbiome/mia?branch=master)

<!-- badges: end -->

This project provides functions and workflows examples for analyses
of microbiome data. The main class for working with microbiome data in this
package is `TreeSummarizedExperiment`.

Currently following things are implemented:

- data wrangling functions (`agglomerate*`, `merge*`, and more)
- CCA analysis via `vegan` package
- Bray-Curtis dissimilarity via `vegan` package
- JSD and UniFrac distance calculation ported from `phyloseq` to work with `TreeSummarizedExperiment` objects
- MDS via the `scater` package for any other distance objects
- import functions for `biom` data, `DADA2` objects, `phyloseq` objects and more

# Contribution

Contributions are welcome in the form of feedback, issues, pull
requests etc, see [contributor guidelines](CONTRIBUTING.md).


## Technical aspects

Let's use a git flow kind of approach. Development version should be done 
against the `master` branch and then merged to `release` for release. 
(https://guides.github.com/introduction/flow/)

## Installation

### Bioc-release

```
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("mia")
```

### Bioc-devel

```
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

# The following initializes usage of Bioc devel
BiocManager::install(version='devel')

BiocManager::install("mia")
```

## How to develop the package

A minimal workflow for package development is as follows.

1. _Fork_ (with Git) mia R package from [https://github.com/microbiome/mia]

1. _Clone_ (with Git) your own fork of mia; or if you have cloned it
  earlier, then make sure that your branch is in sync with the
  _upstream_ repository (see next section)

1. Test package (in R): `devtools::test()` (and fix possible issues)

1. Test examples: `devtools::run_examples()`

1. Update documentation (convert R files to man/ files): `devtools::document()` 

1. Build package: `devtools::build()`

1. Run Bioconductor checks: `BiocCheck::BiocCheck()`

1. Load the updated package: `devtools::load_all()`

1. Install the package: `devtools::install()`

1. Once the package passes all build checks, commit & push to your own
   fork, then make a pull request (PR) to the original
   repository. This will trigger additional automated checks. If they
   fail, inspect the reason, fix accordingly, and update your PR
   (simply commit + push new material to your fork; the PR is
   already open).


Final checks after the R tests; run on command line (replace "R" with
your custom R path if necessary):

- `R CMD build `mia/`
- `R CMD check mia_x.y.z.tar.gz`
- `R CMD BiocCheck mia_x.y.z.tar.gz`
- `R CMD INSTALL mia_x.y.z.tar.gz`


After updating the package, remember to do the following (if applicable):

- update documentation 
- update unit tests (in `tests/`)
- update vignettes (in `vignettes`)
- run all checks from above

After accepted pull request, check if further updates are needed in:

- OMA
- other related packages (e.g. miaViz, miaSim, miaTime)


## Sync with the official Github development version

Before new pull request, check that your own fork is in sync with the
main development version. Command line version goes as follows.

```
# Add the main repository as your _upstream_  
`git remote add upstream git@github.com:microbiome/mia.git`

# Fetch the changes from upstream (main version) and origin (your fork):
git fetch --all

# Merge those changes to your fork:
git merge origin/master
git merge upstream/master

# Finally, commit and push to origin (your version)
# then open a new PR from the Github site if relevant.
git push origin master
```

## Syncing with the official Bioconductor version

This is usually done only by the main developer. In this case, the upstream will be set to Bioconductor: `git remote add upstream git@git.bioconductor.org:packages/mia.git`

- [Workflow to sync between github and bioc versions](https://bioconductor.org/developers/how-to/git/sync-existing-repositories/)

- [Bioconductor Build reports (devel)](https://bioconductor.org/checkResults/devel/bioc-LATEST/)

- [Package page in Bioconductor](https://bioconductor.org/packages/mia/)



# Code of conduct

Please note that the mia project is released with a [Contributor Code of Conduct](https://contributor-covenant.org/version/2/0/CODE_OF_CONDUCT.html).
By contributing to this project, you agree to abide by its terms.

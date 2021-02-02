# mia - Microbiome analysis

<!-- badges: start -->

[![R-CMD-check-bioc-devel](https://github.com/FelixErnst/mia/workflows/R-CMD-check-bioc-devel/badge.svg)](https://github.com/FelixErnst/mia/actions)
[![Codecov test
coverage](https://codecov.io/gh/FelixErnst/mia/branch/master/graph/badge.svg)](https://codecov.io/gh/FelixErnst/mia?branch=master)

<!-- badges: end -->

This project is aimed to provide functions and workflows examples for analyses
of microbiome data. The main class for working with microbiome data in this
package is `TreeSummarizedExperiment`.

Currently following things are implemented:

- data wrangling functions (`agglomerate*`, `merge*`, `meltAssay`)
- JSD and UniFrac distance calculation ported from `phyloseq` to work with `TreeSummarizedExperiment` objects
- CCA analysis via `vegan` package
- Bray-Curtis dissimilarity via `vegan` package
- MDS via the `scater` package
- import functions for `biom` data, `DADA2` objects and `phyloseq` objects

## ToDo

- if `breakaway` ends up on CRAN or Bioconductor, alpha diversity calculation via the `breakaway` package

# Contribution

Feel free to contribute.

## Technical aspects

Let's use a git flow kind of approach. Development version should be done 
against the `master` branch and then merged to `release` for release. 
(https://guides.github.com/introduction/flow/)

# Code of conduct

Please note that the mia project is released with a [Contributor Code of Conduct](https://contributor-covenant.org/version/2/0/CODE_OF_CONDUCT.html).
By contributing to this project, you agree to abide by its terms.

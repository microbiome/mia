# mia - Microbiome analysis

<!-- badges: start -->

[![R-CMD-check-bioc-devel](https://github.com/microbiome/mia/workflows/R-CMD-check-bioc-devel/badge.svg)](https://github.com/microbiome/mia/actions)
[![R-CMD-check-bioc](https://github.com/microbiome/mia/actions/workflows/check-bioc.yml/badge.svg)](https://github.com/microbiome/mia/actions/workflows/check-bioc.yml)
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

Feel free to contribute.

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

# Code of conduct

Please note that the mia project is released with a [Contributor Code of Conduct](https://contributor-covenant.org/version/2/0/CODE_OF_CONDUCT.html).
By contributing to this project, you agree to abide by its terms.

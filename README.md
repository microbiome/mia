# mia - Microbiome analysis <img src="man/figures/mia_logo.png" align="right" width="120" />

<!-- badges: start -->

[![R-CMD-check-bioc-devel](https://github.com/microbiome/mia/workflows/R-CMD-check-bioc-devel/badge.svg)](https://github.com/microbiome/mia/actions)
<!--[![R-CMD-check-bioc](https://github.com/microbiome/mia/actions/workflows/check-bioc.yml/badge.svg)](https://github.com/microbiome/mia/actions/workflows/check-bioc.yml)-->
[![Codecov test
coverage](https://codecov.io/gh/microbiome/mia/branch/master/graph/badge.svg)](https://codecov.io/gh/microbiome/mia?branch=master)

<!-- badges: end -->

## Using the package

This project provides functions and workflows examples for analyses
of microbiome data. The main class for working with microbiome data in this
package is `TreeSummarizedExperiment`. 

For examples of functionality, see the [function reference page](https://microbiome.github.io/mia/reference/index.html).

More information and example workflows are provided in the online
manual [Orchestrating Microbiome Analysis with
Bioconductor](https://microbiome.github.io/OMA).


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


## Contributing

Contributions are welcome in the form of feedback, issues, pull
requests etc, see [contributor guidelines](CONTRIBUTING.md).


## Code of conduct

Please note that the mia project is released with a [Contributor Code of Conduct](https://contributor-covenant.org/version/2/0/CODE_OF_CONDUCT.html).
By contributing to this project, you agree to abide by its terms.

# R package `variantprobs`
[![DOI](https://zenodo.org/badge/184330494.svg)](https://zenodo.org/badge/latestdoi/184330494)

Estimates probabilities and expected numbers of mutations in the tumor genome based on observed mutation frequencies. 

## How to install

### Installing dependencies

`variantprobs` depends on `R` package [edgeR](https://bioconductor.org/packages/release/bioc/html/edgeR.html). To install `edgeR`, run the following commands in `R`:
```{r}
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("edgeR")
```

### Installing `variantprobs`


The easiest way to install `variantprobs` is via `R` package [devtools](https://www.r-project.org/nosvn/pandoc/devtools.html). Run the following commands in `R` to install `devtools`, if it is not already installed:
```{r}
if (!requireNamespace("devtools", quietly = TRUE))
    install.packages("devtools")
```
Then install `variantprobs` as follows:
```{r}
devtools::install_github("c7rishi/variantprobs")
```

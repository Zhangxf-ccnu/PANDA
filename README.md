
<!-- README.md is generated from README.Rmd. Please edit that file -->

# PANDA

<!-- badges: start -->
<!-- badges: end -->

R package supporting the paper “**PANDA enables simultaneous decoding of
cell types and gene expressions for spatial transcriptomics**”.

## Installation

You can install the development version of PANDA like so:

``` r
install.packages("devtools")
devtools::install_github("Zhangxf-ccnu/PANDA")
```

## Usage

This is a basic example which shows you how to use PANDA for
deconvolution.

**Load package**

``` r
library(PANDA)
```

**Load datasets**

-   sc_counts: A matrix of the raw count expression in the scRNA-seq
    reference (cell x gene).
-   sc_labels: A vector of the corresponding cell type labels in the
    scRNA-seq reference.
-   st_counts: A matrix of the raw count expression in the spatial
    transcriptomics data (spot x gene).

**Perform archetypal analysis on the scRNA-seq reference**

``` r
sc_results <- sc_train(sc_counts, sc_labels, n_archetypes_vec = 10)
```

**Perform deconvolution on the spatial transcriptomics data**

``` r
st_results <- st_train(st_counts, sc_results = sc_results)
```

**Extract results**

The cell type proportions for spots can be extracted by

``` r
proportion <- st_results$proportion
```

The cell-type-specific gene expression can be extracted by

``` r
expression <- st_results$mu
```

## Tutorial

A tutorial with examples of the usage of PANDA is available at:

## Contact

Please do not hesitate to contact Mr. Meng-Guo Wang
(<mengguowang@mails.ccnu.edu.cn>) or Dr. Xiao-Fei Zhang
(<zhangxf@ccnu.edu.cn>) to seek any clarifications regarding any
contents or operation of the archive.

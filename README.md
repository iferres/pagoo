<!-- badges: start -->
[![R build status](https://github.com/iferres/pagoo/workflows/R-CMD-check/badge.svg)](https://github.com/iferres/pagoo/actions)
![pkgdown](https://github.com/iferres/pagoo/workflows/pkgdown/badge.svg)
[![codecov](https://codecov.io/gh/iferres/pagoo/branch/master/graph/badge.svg)](https://codecov.io/gh/iferres/pagoo)
<!-- badges: end -->

# Pagoo: An encapsulated OO class system for analyzing bacterial pangenomes in R

*"One object to store them all, one object to find them, one object to query from and with ggplot2 visualize them." (Lord Sauron)*

## Introduction

`pagoo` is an encapsulated, object-oriented class system for analyzing bacterial pangenomes. It uses the [R6](https://r6.r-lib.org/) package as backend. It was designed in order to facilitate and speed-up the comparative analysis of multiple bacterial genomes, standardizing and optimizing routine tasks performed everyday. There are a handful of things done everyday when working with bacterial pangenomes: subset, summarize, extract, visualize and store data. So, `pagoo` is intended to facilitate these tasks as much as possible. 


## Information and Tutorials

Please visit [pagoo webpage](https://iferres.github.io/pagoo/) for more resources.

## Installation

Currently, `pagoo` is only available at GitHub. The easiest way to install this package is by using `devtools`:

``` r
if (!require("devtools")) install.packages("devtools")
devtools::install_github('iferres/pagoo')
```
We have tested `pagoo` in UNIX operating systems. However it should, we cannot warrant it fully works on Windows. Also, in the future we plan to put the code in CRAN and/or Bioconductor repositries.



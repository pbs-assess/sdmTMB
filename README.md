
<!-- README.md is generated from README.Rmd. Please edit that file -->

# sdmTMB

[![Travis build
status](https://travis-ci.org/pbs-assess/sdmTMB.svg?branch=master)](https://travis-ci.org/pbs-assess/sdmTMB)
[![Coverage
status](https://codecov.io/gh/pbs-assess/sdmTMB/branch/master/graph/badge.svg)](https://codecov.io/github/pbs-assess/sdmTMB?branch=master)
[![CRAN
status](https://www.r-pkg.org/badges/version/sdmTMB)](https://cran.r-project.org/package=sdmTMB)
[![Project Status: WIP – Initial development is in progress, but there
has not yet been a stable, usable release suitable for the
public.](https://www.repostatus.org/badges/latest/wip.svg)](https://www.repostatus.org/#wip)

## Installation

You can install sdmTMB with:

``` r
devtools::install_github("pbs-assess/sdmTMB")
```

## Goals of sdmTMB:

  - Be useful for our purpose of fitting high resolution species
    distribution models for groundfish with spatiotemporal GLMMs and
    estimating range shift metrics. Possibly be useful for other
    applications.
  - Have an interface that will be familiar to people who have used
    packages such as glmmTMB, lmer, or mgcv. E.g. accept formulas, use a
    `predict()` and `residuals()` method, and present options as clear
    named arguments.
  - Not fit the kitchen sink, but do what it does elegantly. I.e. be
    opinionated in what it does.

## Example

The main function is `sdmTMB()`. See `?sdmTMB` and `?predict.sdmTMB` for
the most complete examples. More details to come…

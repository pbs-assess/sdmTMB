
<!-- README.md is generated from README.Rmd. Please edit that file -->

# sdmTMB

[![Travis build
status](https://travis-ci.org/pbs-assess/sdmTMB.svg?branch=master)](https://travis-ci.org/pbs-assess/sdmTMB)
[![Project Status: WIP – Initial development is in progress, but there has not yet been a stable, usable release suitable for the public.](https://www.repostatus.org/badges/latest/wip.svg)](https://www.repostatus.org/#wip)

The goal of sdmTMB is to

  - Be useful for our purpose of fitting high resolution species
    distribution models of groundfish with spatiotemporal GLMMs and
    estimating range shift metrics. Possibly be useful for other
    applications.
  - Have an interface that will be familiar to people who have used
    packages such as glmmTMB, lmer, or mgcv. E.g. accept formulas, use a
    `predict()` method, and present options as clear named arguments.
  - Not fit the kitchen sink, but do what it does elegantly. I.e. be
    opinionated in what it does.
  - Maintain clearly documented model code split out into functions as
    much as possible.
  - Adhere to tidyverse package standards where possible.

## Installation

You can install sdmTMB with:

``` r
devtools::install_github("seananderson/sdmTMB")
```

## Example

The main function is `sdmTMB()`. See `?sdmTMB` and `?predict.sdmTMB` for
the most complete examples or the vignette.

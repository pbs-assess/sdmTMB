
<!-- README.md is generated from README.Rmd. Please edit that file -->

# sdmTMB <a href='https://github.com/pbs-assess/sdmTMB'><img src='man/figures/logo-sdmTMB.png' align="right" height="139" /></a>

[![Travis build
status](https://travis-ci.org/pbs-assess/sdmTMB.svg?branch=master)](https://travis-ci.org/pbs-assess/sdmTMB)
[![Project Status: WIP – Initial development is in progress, but there
has not yet been a stable, usable release suitable for the
public.](https://www.repostatus.org/badges/latest/wip.svg)](https://www.repostatus.org/#wip)
[![Coverage
status](https://codecov.io/gh/pbs-assess/sdmTMB/branch/master/graph/badge.svg)](https://codecov.io/github/pbs-assess/sdmTMB?branch=master)

sdmTMB is an R package that implements spatial or spatiotemporal
predictive-process GLMMs (Generalized Linear Mixed Effects Models) using
Template Model Builder ([TMB](https://github.com/kaskr/adcomp)),
[R-INLA](http://www.r-inla.org/), and Gaussian Markov random fields. One
common application is for species distribution models (SDMs).

## Installation

Assuming you have a [C++
compiler](https://support.rstudio.com/hc/en-us/articles/200486498-Package-Development-Prerequisites)
installed, you can install sdmTMB as follows.

First, install INLA:

``` r
install.packages("INLA", repos = c(getOption("repos"), 
  INLA = "https://inla.r-inla-download.org/R/stable"), dep = TRUE)
```

Then install sdmTMB:

``` r
# install.packages("remotes")
remotes::install_github("pbs-assess/sdmTMB")
```

## Functionality

sdmTMB:

  - Fits GLMMs with spatial, spatiotemporal, spatial and spatiotemporal,
    or AR1 spatiotemporal Gaussian Markov random fields with TMB. It can
    also fit spatially varying local trends through time as a random
    field.
  - Uses formula interfaces for fixed effects and any time-varying
    effects (dynamic regression) (e.g. `formula = y ~ 1 + x1,
    time_varying = ~ 0 + x2`), where `y` is the response, `1` represents
    an intercept, `0` omits an intercept, `x1` is a covariate with a
    constant effect, and `x2` is a covariate with a time-varying effect.
  - Can handle formulas with splines from mgcv. E.g., `y ~ s(x, k = 4)`.
  - Uses a `family(link)` format similar to `glm()`, lme4, or glmmTMB.
    This includes Gaussian, Poisson, negative binomial, gamma, binomial,
    lognormal, Student-t, and Tweedie distributions with identity, log,
    inverse, and logit links. E.g., `family = tweedie(link = "log")`.
  - Has `predict()` and `residuals()` methods. The residuals are
    randomized-quantile residuals similar to those implemented in the
    [DHARMa](https://cran.r-project.org/package=DHARMa) package. The
    `predict()` function can take a `newdata` argument similar to `lm()`
    or `glm()` etc. The predictions are bilinear interpolated
    predictive-process predictions (i.e., they make smooth pretty maps).
  - Includes functionality for estimating the centre of gravity or total
    biomass by time step for index standardization.
  - Implements multi-phase estimation for speed.
  - Can optionally allow for anisotropy in the random fields (spatial
    correlation that is directionally dependent).
  - Can generate an SPDE predictive-process mesh based on a clustering
    algorithm and R-INLA or can take any standard R-INLA mesh created
    externally as input.

## Examples

The main function is `sdmTMB()`. See `?sdmTMB` and `?predict.sdmTMB` for
the most complete examples. There is also a simulation function `?sim`
with some examples. There are some vignettes you can see if you build
with `devtools::install_github("pbs-assess/sdmTMB", build_vignettes =
TRUE)` or look
[here](https://github.com/pbs-assess/sdmTMB/tree/master/vignettes).

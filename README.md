
<!-- README.md is generated from README.Rmd. Please edit that file -->

# sdmTMB <a href='https://github.com/pbs-assess/sdmTMB'><img src='man/figures/logo-sdmTMB.png' align="right" height="139" /></a>

> Spatial and spatiotemporal GLMMs with TMB

<!-- badges: start -->
[![R-CMD-check](https://github.com/pbs-assess/sdmTMB/workflows/R-CMD-check/badge.svg)](https://github.com/pbs-assess/sdmTMB/actions)
[![Project Status: WIP – Initial development is in progress, but there has not yet been a stable, usable release suitable for the public.](https://www.repostatus.org/badges/latest/wip.svg)](https://www.repostatus.org/#wip)
[![Coverage
status](https://codecov.io/gh/pbs-assess/sdmTMB/branch/master/graph/badge.svg)](https://codecov.io/github/pbs-assess/sdmTMB?branch=master)
<!-- badges: end -->

sdmTMB is an R package that fits spatial and spatiotemporal GLMMs (Generalized Linear Mixed Effects Models) using Template Model Builder ([TMB](https://github.com/kaskr/adcomp)), [R-INLA](https://www.r-inla.org/), and Gaussian Markov random fields. One common application is for species distribution models (SDMs).

## Table of contents

-   [Installation](#installation)
-   [Functionality](#functionality)
-   [Basic use](#basic-use)
-   [Advanced functionality](#advanced-functionality)
    -   [Time-varying coefficients](#time-varying-coefficients)
    -   [Spatially varying coefficients
        (SVC)](#spatially-varying-coefficients-svc)
    -   [Simulating data](#simulating-data)
    -   [Bayesian MCMC sampling with
        Stan](#bayesian-mcmc-sampling-with-stan)
    -   [Sampling from the joint precision
        matrix](#sampling-from-the-joint-precision-matrix)
    -   [Calculating uncertainty on spatial
        predictions](#calculating-uncertainty-on-spatial-predictions)
    -   [Priors](#priors)
    -   [Cross validation](#cross-validation)
    -   [Turning off random fields](#turning-off-random-fields)
    -   [Using a custom INLA mesh](#using-a-custom-inla-mesh)

## Installation

Assuming you have a [C++
compiler](https://support.rstudio.com/hc/en-us/articles/200486498-Package-Development-Prerequisites)
installed, you can install sdmTMB:

``` r
# install.packages("remotes")
remotes::install_github("pbs-assess/sdmTMB")
```

## Functionality

sdmTMB:

-   Fits GLMMs with spatial, spatiotemporal, spatial and spatiotemporal,
    or AR1 spatiotemporal Gaussian Markov random fields with TMB.
-   Uses formula interfaces for fixed effects and any time-varying
    effects (dynamic regression)
    (e.g. `formula = y ~ 1 + x1 + (1 | g), time_varying = ~ 0 + x2`),
    where `y` is the response, `1` represents an intercept, `0` omits an
    intercept, `x1` is a covariate with a constant effect, `(1 | g)` is
    a random intercept across groups `g`, and `x2` is a covariate with a
    time-varying effect.
-   Can fit spatially varying coefficients as a random field
    (e.g. `spatial_varying = ~ 0 + x3`).
-   Can handle GAMs (generalized additive models) with penalized
    smoothers from mgcv. E.g., `y ~ s(x)`.
-   Can handle linear breakpoint or logistic threshold fixed effects:
    `y ~ breakpt(x1)` or `y ~ logistic(x2)`.
-   Uses a `family(link)` format similar to `glm()`, lme4, or glmmTMB.
    This includes Gaussian, Poisson, negative binomial, gamma, binomial,
    lognormal, Student-t, and Tweedie distributions with identity, log,
    inverse, and logit links. E.g., `family = tweedie(link = "log")`.
-   Has `predict()` and `residuals()` methods. The residuals are
    randomized-quantile residuals similar to those implemented in the
    [DHARMa](https://cran.r-project.org/package=DHARMa) package. The
    `predict()` function can take a `newdata` argument similar to `lm()`
    or `glm()` etc. The predictions are bilinear interpolated
    predictive-process predictions (i.e., they make smooth pretty maps).
-   Has a simulation function `simulate()` for simulating from existing
    fits (e.g., for DHARMa), `sdmTMB_simulate()` for generating
    simulated data from scratch, and `sdmTMB_cv()` for cross-validation
    testing of model accuracy or comparing across model configurations.
-   Includes functionality for estimating the centre of gravity or total
    biomass by time step for index standardization.
-   Can optionally allow for anisotropy in the random fields (spatial
    correlation that is directionally dependent) and barriers (e.g.,
    land for ocean species) to spatial correlation.
-   Can interpolate over missing time slices or forecast onto future
    time slices.
-   Can generate an SPDE predictive-process mesh or can take any
    standard R-INLA mesh created externally as input.

See `?sdmTMB` and `?predict.sdmTMB` for the most complete examples. Also
see the vignettes (‘Articles’) on the [documentation
site](https://pbs-assess.github.io/sdmTMB/index.html).

## Basic use

An sdmTMB model requires a data frame that contains a response column,
columns for any predictors, and columns for spatial coordinates. It
usually makes sense to convert the spatial coordinates to an equidistant
projection such as UTMs such that distance remains constant throughout
the study region \[e.g., using `sf::st_transform()`\]. Here, we
illustrate a spatial model fit to Pacific cod (*Gadus macrocephalus*)
trawl survey data from Queen Charlotte Sound, BC, Canada. Our model
contains a main effect of depth as a penalized smoother, a spatial
random field, and Tweedie observation error. Our data frame `pcod`
(built into the package) has a column `year` for the year of the survey,
`density` for density of Pacific cod in a given survey tow, `present`
for whether `density > 0`,
depth`for depth in meters of that tow, and spatial coordinates`X`and`Y\`,
which are UTM coordinates in kilometres.

``` r
library(dplyr)
library(ggplot2)
library(sdmTMB)
head(pcod)
```

    #> # A tibble: 3 × 6
    #>    year density present depth     X     Y
    #>   <int>   <dbl>   <dbl> <dbl> <dbl> <dbl>
    #> 1  2003   113.        1   201  446. 5793.
    #> 2  2003    41.7       1   212  446. 5800.
    #> 3  2003     0         0   220  449. 5802.

We start by creating a mesh object that contains matrices to apply the
SPDE approach.

``` r
mesh <- make_mesh(pcod, xy_cols = c("X", "Y"), cutoff = 10)
```

Here, `cutoff` defines the minimum allowed distance between points in
the units of `X` and `Y` (km). Alternatively, we could have created any
mesh via the INLA package and supplied it to `make_mesh()`. We can
inspect our mesh object with the associated plotting method
`plot(mesh)`.

Fit a spatial model with a smoother for depth:

``` r
fit <- sdmTMB(
  density ~ s(depth, k = 5),
  data = pcod,
  family = tweedie(link = "log"),
  mesh = mesh,
  spatial = "on"
)
```

Print the model fit:

``` r
fit
#> Spatiotemporal model fit by ML ['sdmTMB']
#> Formula: density ~ s(depth, k = 5)
#> Mesh: mesh
#> Data: pcod
#> Family: tweedie(link = 'log')
#>             coef.est coef.se
#> (Intercept)     2.34    0.23
#> sdepth        -26.96   21.52
#> 
#> Smooth terms:
#>            Std. Dev.
#> sds(depth)     34.03
#> 
#> Dispersion parameter: 12.70
#> Tweedie p: 1.58
#> Matern range: 17.47
#> Spatial SD: 1.86
#> ML criterion at convergence: 6403.637
#> 
#> See ?tidy.sdmTMB to extract these values as a data frame.
```

The output indicates our model was fit by maximum (marginal) likelihood
(`ML`). We also see the formula, mesh, fitted data, and family. Next we
see any estimated main effects including the linear component of the
smoother (`sdepth`), the standard deviation on the smoother weights
(`sds(depth)`), the Tweedie dispersion and power parameters, the Matérn
range distance (distance at which points are effectively indepenent),
the marginal spatial field standard deviation, and the negative log
likelihood at convergence.

We can extract parameters as a data frame:

``` r
tidy(fit, conf.int = TRUE)
#> # A tibble: 1 × 5
#>   term        estimate std.error conf.low conf.high
#>   <chr>          <dbl>     <dbl>    <dbl>     <dbl>
#> 1 (Intercept)     2.34     0.227     1.89      2.78
tidy(fit, effects = "ran_pars", conf.int = TRUE)
#> # A tibble: 4 × 5
#>   term      estimate std.error conf.low conf.high
#>   <chr>        <dbl> <lgl>        <dbl>     <dbl>
#> 1 range        17.5  NA           10.5      29.1 
#> 2 phi          12.7  NA           11.9      13.5 
#> 3 sigma_O       1.86 NA            1.51      2.29
#> 4 tweedie_p     1.58 NA            1.57      1.60
```

Plot the smoother effect:

``` r
plot_smooth(fit, ggplot = TRUE)
```

<img src="man/figures/README-unnamed-chunk-8-1.png" width="50%" />

Predict on new data:

``` r
p <- predict(fit, newdata = qcs_grid)
```

``` r
head(p)
```

    #> # A tibble: 3 × 7
    #>       X     Y depth   est est_non_rf est_rf omega_s
    #>   <dbl> <dbl> <dbl> <dbl>      <dbl>  <dbl>   <dbl>
    #> 1   456  5636  347. -2.32      -2.34 0.0155  0.0155
    #> 2   458  5636  223.  2.06       2.02 0.0475  0.0475
    #> 3   460  5636  204.  3.10       3.02 0.0796  0.0796

``` r
ggplot(p, aes(X, Y, fill = exp(est))) + geom_raster() +
  scale_fill_viridis_c(trans = "sqrt")
```

<img src="man/figures/README-unnamed-chunk-12-1.png" width="50%" />

We could switch to a presence-absence model by changing the response
column and family:

``` r
fit <- sdmTMB(
  present ~ s(depth, k = 5),
  data = pcod,
  family = binomial(link = "logit"),
  mesh = mesh
)
```

We could instead fit a spatiotemporal model by specifying the `time`
column and a spatiotemporal structure:

``` r
fit_spatiotemporal <- sdmTMB(
  density ~ s(depth, k = 5), 
  family = tweedie(link = "log"), 
  data = pcod, 
  mesh = mesh,
  time = "year", 
  spatial = "off", 
  spatiotemporal = "ar1"
)
```

If we wanted to create an area-weighted standardized population index,
we could predict on a grid covering the entire survey (`qcs_grid`) with
grid cell area 4 (2 x 2 km) and pass the predictions to `get_index()`:

``` r
p_st <- predict(fit_spatiotemporal, newdata = qcs_grid, 
  return_tmb_object = TRUE, area = 4)
index <- get_index(p_st)
ggplot(index, aes(year, est)) +
  geom_ribbon(aes(ymin = lwr, ymax = upr), fill = "grey90") +
  geom_line(lwd = 1, colour = "grey30") +
  labs(x = "Year", y = "Biomass (kg)")
```

<img src="man/figures/README-unnamed-chunk-15-1.png" width="50%" />

Or the center of gravity:

``` r
cog <- get_cog(p_st, format = "wide")
ggplot(cog, aes(est_x, est_y, colour = year)) +
  geom_pointrange(aes(xmin = lwr_x, xmax = upr_x)) +
  geom_pointrange(aes(ymin = lwr_y, ymax = upr_y)) +
  scale_colour_viridis_c()
```

<img src="man/figures/README-unnamed-chunk-16-1.png" width="50%" />

## Advanced functionality

### Time-varying coefficients

Time-varying (random walk) effect of depth:

``` r
mesh <- make_mesh(pcod, c("X", "Y"), cutoff = 10)
fit <- sdmTMB(
  density ~ 1, 
  time_varying = ~ 0 + depth,
  family = tweedie(link = "log"),
  data = pcod,
  mesh = mesh,
  spatial = "off",
  spatiotemporal = "ar1",
  time = "year",
  silent = FALSE # see progress
)
```

Time-varying intercept:

``` r
fit <- sdmTMB(
  density ~ 0 + s(depth, k = 5), 
  time_varying = ~ 1, 
  data = pcod,
  mesh = mesh,
  time = "year",
  silent = FALSE
)
```

### Spatially varying coefficients (SVC)

Spatially varying effect of time:

``` r
pcod$year_scaled <- as.numeric(scale(pcod$year))
fit <- sdmTMB(
  density ~ s(depth, k = 5) + year_scaled,
  data = pcod,
  mesh = mesh, 
  family = tweedie(link = "log"),
  spatial_varying = ~ 0 + year_scaled, 
  time = "year",
  spatiotemporal = "off"
)
```

See `zeta_s` in the output, which represents the coefficient varying in
space. You’ll want to ensure you set up your model such that it ballpark
has a mean of 0 (e.g., by including it in `formula` too).

``` r
qcs_grid$year_scaled <- (qcs_grid$year - mean(pcod$year)) / sd(pcod$year)
p <- predict(fit, newdata = subset(qcs_grid, year == 2011))
ggplot(p, aes(X, Y, fill = zeta_s)) + geom_raster() +
  scale_fill_gradient2()
```

<img src="man/figures/README-unnamed-chunk-20-1.png" width="50%" />

### Simulating data

#### Simulating data from scratch:

``` r
predictor_dat <- expand.grid(
  X = seq(0, 1, length.out = 100), Y = seq(0, 1, length.out = 100)
)
mesh <- make_mesh(predictor_dat, xy_cols = c("X", "Y"), cutoff = 0.05)
sim_dat <- sdmTMB_simulate(
  formula = ~ 1,
  data = predictor_dat,
  mesh = mesh,
  family = poisson(link = "log"),
  range = 0.3,
  sigma_O = 0.4,
  seed = 1,
  B = 1 # B0 = intercept
)
head(sim_dat)
#> # A tibble: 6 × 7
#>        X     Y omega_s    mu   eta observed `(Intercept)`
#>    <dbl> <dbl>   <dbl> <dbl> <dbl>    <dbl>         <dbl>
#> 1 0          0  -0.154  2.33 0.846        1             1
#> 2 0.0101     0  -0.197  2.23 0.803        0             1
#> 3 0.0202     0  -0.240  2.14 0.760        2             1
#> 4 0.0303     0  -0.282  2.05 0.718        2             1
#> 5 0.0404     0  -0.325  1.96 0.675        3             1
#> 6 0.0505     0  -0.367  1.88 0.633        2             1

# sample 200 points for fitting:
set.seed(1)
sim_dat_obs <- sim_dat[sample(seq_len(nrow(sim_dat)), 200), ]
```

``` r
ggplot(sim_dat, aes(X, Y)) +
  geom_raster(aes(fill = exp(eta))) + # mean without observation error
  geom_point(aes(size = observed), data = sim_dat_obs, pch = 21) +
  scale_fill_viridis_c() +
  scale_size_area() +
  coord_cartesian(expand = FALSE)
```

<img src="man/figures/README-unnamed-chunk-22-1.png" width="50%" />

Fit to the simulated data:

``` r
mesh <- make_mesh(sim_dat_obs, xy_cols = c("X", "Y"), cutoff = 0.05)
fit <- sdmTMB(
  observed ~ 1,
  data = sim_dat_obs,
  family = poisson(),
  mesh = mesh
)
```

#### Simulating from an existing fit

``` r
s <- simulate(fit, nsim = 100)
dim(s)
#> [1] 200 100
s[1:3,1:4]
#>      [,1] [,2] [,3] [,4]
#> [1,]    5    2    4    4
#> [2,]    2    2    3    2
#> [3,]    5    5    2    4
```

Using those simulations to check DHARMa residuals:

``` r
res <- DHARMa::createDHARMa(
  simulatedResponse = s,
  observedResponse = sim_dat_obs$observed,
  fittedPredictedResponse = predict(fit)$est_non_rf # identity link
)
DHARMa::plotQQunif(res, testUniformity = FALSE, testOutliers = FALSE, testDispersion = FALSE)
```

<img src="man/figures/README-unnamed-chunk-25-1.png" width="50%" />

Comparing those to the built-in randomized-quantile residuals:

``` r
r <- residuals(fit)
qqnorm(r);qqline(r)
```

<img src="man/figures/README-unnamed-chunk-26-1.png" width="50%" />

### Bayesian MCMC sampling with Stan

The fitted model can be passed to the tmbstan package to sample from the
posterior with Stan. Note this can be slow for larger models.

``` r
# only 1 chain and 400 iterations for speed:
fit_mcmc <- tmbstan::tmbstan(fit$tmb_obj, chains = 1, iter = 400)
```

Internal paramater posteriors:

``` r
print(fit_mcmc, pars = c("b_j", "omega_s[1]"))
#> Inference for Stan model: sdmTMB.
#> 1 chains, each with iter=400; warmup=200; thin=1; 
#> post-warmup draws per chain=200, total post-warmup draws=200.
#> 
#>             mean se_mean   sd  2.5%   25%   50%  75% 97.5% n_eff Rhat
#> b_j         0.99    0.03 0.15  0.62  0.93  1.00 1.06  1.27    35 1.00
#> omega_s[1] -0.07    0.03 0.23 -0.50 -0.23 -0.06 0.10  0.33    63 1.01
#> 
#> Samples were drawn using NUTS(diag_e) at Mon Dec 20 16:19:18 2021.
#> For each parameter, n_eff is a crude measure of effective sample size,
#> and Rhat is the potential scale reduction factor on split chains (at 
#> convergence, Rhat=1).
```

Predicting with the Stan/tmbstan model:

``` r
pred_mcmc <- predict(fit, newdata = qcs_grid, tmbstan_model = fit_mcmc)
# Each row has 200 posterior samples for a row of the `newdata` data frame:
dim(pred_mcmc)
#> [1] 65826   200
```

### Sampling from the joint precision matrix

We can take samples from the implied parameter distribution assuming an
MVN covariance matrix on the internal parameterization:

``` r
samps <- gather_sims(fit, n_sims = 1000)
ggplot(samps, aes(.value)) + geom_histogram() +
  facet_wrap(~.variable, scales = "free_x")
#> `stat_bin()` using `bins = 30`. Pick better value with `binwidth`.
```

<img src="man/figures/README-unnamed-chunk-30-1.png" width="50%" />

### Calculating uncertainty on spatial predictions

The fastest way to get pointwise prediction uncertainty is to use the
MVN samples:

``` r
p <- predict(fit, newdata = predictor_dat, sims = 500)
predictor_dat$se <- apply(p, 1, sd)
ggplot(predictor_dat, aes(X, Y, fill = se)) +
  geom_raster() +
  scale_fill_viridis_c(option = "A") +
  coord_cartesian(expand = FALSE)
```

<img src="man/figures/README-unnamed-chunk-31-1.png" width="50%" />

### Priors

Priors (or ‘penalties’) can be placed on most parameters. For example,
here we place a PC prior on the Matérn random field parameters, a
standard normal prior on the effect of depth, a Normal(0, 10^2) prior on
the intercept, and a half-normal prior on the Tweedie dispersion
parameter (`phi`):

``` r
mesh <- make_mesh(pcod, c("X", "Y"), cutoff = 10)
fit <- sdmTMB(
  density ~ depth_scaled,
  data = pcod,
  family = tweedie(),
  mesh = mesh,
  priors = sdmTMBpriors(
    matern_s = pc_matern(range_gt = 10, sigma_lt = 5),
    b = normal(c(0, 0), c(1, 10)),
    phi = halfnormal(0, 15)
  )
)
```

### Cross validation

sdmTMB has built-in functionality for cross-validation. If we were to
set a `future::plan()`, the folds would be fit in parallel:

``` r
mesh <- make_mesh(pcod, c("X", "Y"), cutoff = 10)
## Set parallel processing if desired:
# library(future)
# plan(multisession)
m_cv <- sdmTMB_cv(
  density ~ s(depth, k = 5),
  data = pcod, mesh = mesh,
  family = tweedie(link = "log"), k_folds = 2
)
#> Running fits with `future.apply()`.
#> Set a parallel `future::plan()` to use parallel processing.
# Sum of log likelihoods of left-out data:
m_cv$sum_loglik
#> [1] -6806.148
# Expected log pointwise predictive density from left-out data:
# (average likelihood density)
m_cv$elpd
#> [1] -1.009422
```

### Turning off random fields

We can turn off the random fields for model comparison:

``` r
fit_sdmTMB <- sdmTMB(
  present ~ poly(depth_scaled, 2),
  data = pcod, mesh = mesh,
  spatial = "off",
  family = binomial()
)
fit_glm <- glm(
  present ~ poly(depth_scaled, 2),
  data = pcod,
  family = binomial()
)

tidy(fit_sdmTMB)
#> # A tibble: 3 × 3
#>   term                   estimate std.error
#>   <chr>                     <dbl>     <dbl>
#> 1 (Intercept)              -0.426    0.0573
#> 2 poly(depth_scaled, 2)1  -31.7      3.03  
#> 3 poly(depth_scaled, 2)2  -66.9      4.09
broom::tidy(fit_glm)
#> # A tibble: 3 × 5
#>   term                   estimate std.error statistic  p.value
#>   <chr>                     <dbl>     <dbl>     <dbl>    <dbl>
#> 1 (Intercept)              -0.426    0.0573     -7.44 1.03e-13
#> 2 poly(depth_scaled, 2)1  -31.7      3.03      -10.5  1.20e-25
#> 3 poly(depth_scaled, 2)2  -66.9      4.09      -16.4  3.50e-60
```

### Using a custom INLA mesh

Defining a mesh directly with INLA:

``` r
bnd <- INLA::inla.nonconvex.hull(cbind(pcod$X, pcod$Y), convex = -0.1)
mesh_inla <- INLA::inla.mesh.2d(
  boundary = bnd,
  max.edge = c(25, 50)
)
mesh <- make_mesh(pcod, c("X", "Y"), mesh = mesh_inla)
plot(mesh)
```

<img src="man/figures/README-unnamed-chunk-35-1.png" width="50%" />

``` r
fit <- sdmTMB(
  density ~ s(depth, k = 5),
  data = pcod,
  family = tweedie(link = "log"),
  mesh = mesh
)
```

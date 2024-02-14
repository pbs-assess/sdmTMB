# sdmTMB (development version)

* Issue error if `time` column has NAs. #298 #299

* Fix bug in `get_cog(..., format = "wide")` where the time column was
  hardcoded to `"year"` by accident.

* Poisson-link delta models now use a `type` argument in `delta_gamma()` and
  `delta_lognormal()`. `delta_poisson_link_gamma()` and
  `delta_poisson_link_lognormal()` are deprecated. #290
  
* Delta families can now pass links that are different from the default 
  `"logit"` and `"log"`. #290

# sdmTMB 0.4.2

* Force rebuild of CRAN binaries to fix issue with breaking Matrix ABI change
  causing `NaN gradient` errors. #288 #287

* Fix crash in if `sdmTMB(..., do_index = TRUE)` and `extra_time` supplied along
  with `predict_args = list(newdata = ...)` that lacked `extra_time` elements.

* Allow `get_index()` to work with missing time elements.

* Add the ability to pass a custom randomized quantile function `qres_func`
  to `residuals.sdmTMB()`.

* Add check for factor random intercept columns in `newdata` to avoid a crash.
  #278 #280

* Improve warnings/errors around use of `do_index = TRUE` and `get_index()`
  if `newdata = NULL`. #276

* Fix prediction with `offset` when `newdata` is `NULL` but `offset` is
  specified. #274

* Fix prediction failure when both `offset` and `nsim` are provided and
  model includes `extra_time`. #273

# sdmTMB 0.4.1

* Fix memory issues detected by CRAN 'Additional issues' clang-UBSAN, valgrind.

* Fix a bug predicting on new data with a specified offset and `extra_time`. 
  #270

* Add warning around non-factor handling of the `spatial_varying` formula. #269

* Add experimental `set_delta_model()` for plotting delta models with
  `ggeffects::ggpredict()` (GitHub version only until next CRAN version).

# sdmTMB 0.4.0

* Move add_barrier_mesh() to sdmTMBextra to avoid final INLA dependency.
  https://github.com/pbs-assess/sdmTMBextra

* Switch to using the new fmesher package for all mesh/SPDE calculations. INLA
  is no longer a dependency.
  
* Switch to `diagonal.penalty = FALSE` in mgcv::smoothCon(). 
  This changes the scale of the linear component of the smoother, but
  should result in the same model.
  https://github.com/glmmTMB/glmmTMB/issues/928#issuecomment-1642862066
  
* Implement cross validation for delta models #239

* Remove ELPD from cross validation output. Use sum_loglik instead. #235

* Turn on Newton optimization by default. #182

* print() now checks sanity() and issues a warning if there may be issues. #176

* Poisson-link delta models and censored likelihood distributions have been made
  considerably more robust. #186
  
* Standard errors are now available on SD parameters etc. in tidy() #240

* Fix bug in print()/tidy() for delta-model positive model component sigma_E.
  A recently introduce bug was causing sigma_E for the 2nd model to be reported
  as the 1st model component sigma_E.
  
* Add new anisotropy plotting function.

* Add anisotropic range printing. #149 by @jdunic

# sdmTMB 0.3.0

* Create the sdmTMBextra package to remove rstan/tmbstan helpers, which
  were causing memory sanitizer errors on CRAN.
  https://github.com/pbs-assess/sdmTMBextra
  
* The following functions are affected:

  - `predict.sdmTMB()` now takes `mcmc_samples`, which is output from
    `sdmTMBextra::extract_mcmc()`.
  - `simulate.sdmTMB()` now takes `mcmc_samples`, which is output from
    `sdmTMBextra::extract_mcmc()`.
  - `residuals.sdmTMB()` now takes `mcmc_samples`, which is output
    `sdmTMBextra::predict_mle_mcmc()`. This only affects 
    `residuals(..., type = "mle-mcmc")`.

* Move `dharma_residuals()` to 
  [sdmTMBextra](https://github.com/pbs-assess/sdmTMBextra) to reduce heavy
  dependencies.
  
* See examples in the Bayesian and residuals vignettes or in the help files for
  those functions within sdmTMBextra.

# sdmTMB 0.2.2

* Various fixes to pass CRAN checks. #158

* Fix memory issue highlighted by Additional issues CRAN checks. #158

* 'offset' argument can now be a character value indicating a column name. This
  is the preferred way of using an offset with parallel cross validation. #165

* Fix parallel cross validation when using an offset vector. #165

* Add leave-future-out cross validation functionality. #156

* Example data `qcs_grid` is no longer replicated by year to save package
  space. #158

* Add message with `tidy(fit, "ran_pars")` about why SEs are NA.

* Add anisotropy to `print()` #157

* Fix `predict(..., type = "response", se_fit = TRUE)`, which involves issuing
  a warning and sticking to link space. #140
  
# sdmTMB 0.2.1

* Fixes for resubmission to CRAN.

# sdmTMB 0.2.0

* Initial submission to CRAN.

# sdmTMB 0.1.4

* Relax range parameter `sanity()` check from 1x to 1.5x the greatest
  distance in the data.

* Add Pearson residuals for several families. `residuals(fit, type = "pearson")`
  Useful for checking for overdispersion with N > 1 binomial or Poisson
  families, among other uses. See the `overdisp_fun()` function at:
  https://bbolker.github.io/mixedmodels-misc/glmmFAQ.html#testing-for-overdispersioncomputing-overdispersion-factor

* Fix bug when using `residuals()` or `simulate()` with binomial families
  specified via `cbind()` or `weights = N`. The binomial sample size wasn't
  being passed through typically resulting in Inf/-Inf.

* Add mixture families: `gamma_mix()`, `lognormal_mix()` and associated
  delta/hurdle families: `delta_gamma_mix()`, `delta_lognormal_mix()`. These
  families feature a mixture of two distributions with different means but
  shared variance parameters.
  
* Add `delta_beta()` family.

# sdmTMB 0.1.3

* Tweak `sanity()` checking of standard error size.

* Export previously experimental `plot_anisotropy()` function. The old function
is now `plot_anisotropy2()`.

* Allow passing offset data through `predict.sdmTMB()` via `offset` argument.

# sdmTMB 0.1.2

* Switch `effects = 'ran_vals'` for random intercept values from `tidy.sdmTMB()`
  to match the broom.mixed package.

* Make `tidy.sdmTMB()` return a tibble if the tibble package is installed. Note 
  this could affect old code since `drop = FALSE` is the default for tibbles
  but `drop = TRUE` is the default for data frames (i.e., tibbles always return
  a data frame when subsetted).

* Fix longstanding issue with predicting on newdata with mgcv's `t2()`. 
  Previously this was disabled because of issues. It now works as expected.
  
* Add `knots` argument in `sdmTMB()`, which is passed to mgcv. A common use
  would be to specify end points in a cyclical spline 
  (e.g., `s(x, bs = 'cc', k = 4), knots = list(x = c(1, 3, 5, 7))`) when the
  data don't extend fully to the boundaries that should match up.

# sdmTMB 0.1.1

* Preparing for release on CRAN.

* Add time-varying AR1 option (originally was always a random walk). See 
  `time_varying_type` argument in `?sdmTMB`.

* Allow prediction on `newdata` with missing time elements. #130

* Add check for `offset()` (which *does not* work in sdmTMB, use the `offset`
  argument instead). #131
  
* Add check for random slopes (sdmTMB currently only does random intercepts,
  although slopes can vary spatially). #131

# sdmTMB 0.1.0

* ADREPORT several parameters in natural space.
  <https://github.com/pbs-assess/sdmTMB/discussions/113>

* Improve robustness of model `print()` to more esoteric mgcv smoothers.

* Let `sims_var` work with multiple spatially varying slopes (`zeta_s`); return
  output in named list by coefficients. #107
  
* Add `threshold_coefs` to `sdmTMB_simulate()`.

* Don't make a fake mesh for non-spatial model (faster).

# sdmTMB 0.0.26.9001

* Add vignettes on visreg, ggeffects, and delta families (thanks J. Indivero!)
  #83 #87 #89 Forecasting and presence-only vignettes to be merged in soon.

* Add support for emmeans package. See `?emmeans.sdmTMB` for examples.

* Add support for effects package. The `ggeffects::ggeffect()` function
  can be used to make fast marginal effects plots. `ggeffects::ggpredict()`
  works with a custom fork of ggeffects. A pull request will be made shortly.
  #101 

* Add `vcov()`, `fixef()`, `df.residual`(), `formula()`, `terms()`, and
  `model.frame()` methods.

* Add support for `"cloglog"` link. Code adapted from glmmTMB for robust
  likelihood implementation.

* For delta models, by default share the anisotropy parameters as in VAST. 
  Separate anisotropy (old behavior) can be estimated with
  `control = sdmTMBcontrol(map = list(ln_H_input = factor(c(1, 2, 3, 4))))`
  
* Add experimental `do_index`, `predict_args`, and `index_args` in `sdmTMB()`.
  These can be used to perform prediction and index calculation at the same
  time as fitting. For very large datasets or meshes this can save time
  compared to fitting, predicting, and index calculation in 3 separate steps
  since the TMB AD object doesn't have to be rebuilt. This will somewhat slow
  down the initial fitting.
  
* Remove `max_gradient` and `bad_eig` from `get_index()` output.

* Use unique locations on prediction for huge speedups on large `newdata`
  gridded data.

* Fix bug where in rare cases `get_index()` would return gibberish small values.

* Add `bayesian` argument, which when `TRUE` adds Jacobian adjustments for
  non-linear transformed parameters. This should be `TRUE` if the model 
  will be passed to tmbstan, but `FALSE` otherwise. #95
  
* Add experimental and not-yet-exported `sdmTMB:::plot_anisotropy2()`.

* Add many anisotropy, delta model, and index calculation unit tests.

# sdmTMB 0.0.24.9001

* Enable random walk random field TMB simulation in `sdmTMB_simulate()`.

* Add check for irregular time with AR1 or random walk processes.

* Fix bugs introduced by delta model code (offsets with `extra_time` and
  threshold model prediction).
  
* Fix bug in `sanity()` message with small random field SDs.

# sdmTMB 0.0.24.9000

* Add support for 'delta' (or 'hurdle') models. See examples and documentation
  in `?sdmTMB`. This has resulted in a substantial restructuring of the
  internal model code. By default both model components (e.g., binomial & Gamma)
  share the same formula, spatial, and spatiotemporal structure, but these
  can be separated by supplying argument values in lists where the first
  element corresponds to the first model and the second element corresponds to
  the second model (with some limitations as described in `?sdmTMB`
  documentation 'Details').

* Add support for multiple spatially varying coefficients (used to be limited to
  a single variable).

* Add compatibility with the 'visreg' package for visualizing conditional
  effects of parameters. See `?visreg_delta` for examples.

* Add MCMC residual type to `residuals.sdmTMB()`. These are a 'better' residuals
  but slower to calculate. See documentation 'Details' in `?residuals.sdmTMB`.

* Make `offset` an argument in `sdmTMB()`. Using the reserved word `offset` in
  the formula is now deprecated.

* Add `sanity()` function to perform some basic sanity checks on model fits.

* Make an `sdmTMB()` model object compatible with `update()` method.

* Remove several deprecated arguments.

* Overhaul examples in `?sdmTMB`.

* Use faster "low-rank sparse hessian bias-correction" TMB bias correction.

* Add parallel processing support. See `parallel` argument in `sdmTMBcontrol`.
  By default, grabs value of `sdmTMB.cores` option. E.g.
  `options(sdmTMB.cores = 4)`. Only currently enabled on Mac/Linux.
  Using too many cores can be much slower than 1 core.
  
* Use 'cli' package `cli_abort()`/`cli_warn()`/`cli_inform()` over
  `stop()`/`warning()`/`message()`.

* Add many unit tests.

# sdmTMB 0.0.23.9000

* A package version number that was used for internal testing in the 'delta'
  branch by several people.

# sdmTMB 0.0.22.9001

* Switch to TMBad library for ~3-fold speedup(!)

# sdmTMB 0.0.22.9000

* Fix bug in predictions with `poly(..., raw = FALSE)` on newdata. #77

# sdmTMB 0.0.21.9009

* Add experimental `sdmTMB_stacking()` for ensemble model stacking weights.

* Add fake mesh if random fields are all off. #59

* Make `predict(..., newdata = NULL)` also use `last.par.best` instead of
  `last.par` to match `newdata = df`.

* Fix bug in MVN fixed-effect prior indexing

* `sims` and `n_sims` arguments have been deprecated and standardized
  to `nsim` to match the `simulate()` S3 method.

* Bias correction on `get_index()` and `get_cog()` is now selective and 
  is just applied to the necessary derived parameters.

* INLA projection matrix 'A' is now shared across spatial and spatiotemporal
  fields.

* Add `add_utm_columns()` to ease adding UTM columns.

# sdmTMB 0.0.20.9001

* Add `dharma_residuals()`.

* Fix bug in `simulate.sdmTMB()` and `residuals.sdmTMB()` for binomial family.

# sdmTMB 0.0.20.9000

*  Smoothers now appear in `print()` output. The format should roughly match brms.
   The main-effect component (e.g., `sdepth` for `s(depth)`) represents the
   linear component and the random effect (e.g., `sds(depth)`) component in
   the output corresponds to the standard deviation of the penalized weights.

*  Add `censored_poisson(link = 'log')` family; implemented by @joenomiddlename

* `fields` in `sdmTMB()` is now deprecated and replaced by `spatiotemporal`.

* `include_spatial` in `sdmTMB()` is now deprecated and replaced by `spatial`.

* `spatial_only` in `sdmTMB()` is now deprecated and replaced by `spatiotemporal`. 
   E.g. `spatial_only = TRUE` is now `spatiotemporal = 'off'` or leaving 
   `time = NULL`.
   
* `spde` in `sdmTMB()` is now deprecated and replaced by `mesh`.

* `sdmTMB_simulate()` is new and will likely replace `sdmTMB_sim()` eventually. 
  `sdmTMB_simulate()` is set up to take a formula and a data frame and is easier
  to use if you want different spatial observations (and covariates) for each
  time slice. It can also take a fitted model and modify parts of it to simulate.
  Finally, this function uses TMB for simulation and so is much faster and
  more flexible in what it can simulate (e.g., anisotropy) than the previous version.
  
* `spatial_trend` is now `spatial_varying` and accepts a one-sided formula
  *with a single predictor* of any coefficient that should varying in space as a
  random field. Note that you may want to include a fixed effect for the same
  variable to improve interpretability. If the (scaled) time column is used, it will
  represent a local-time-trend model as before.
  
* The Tweedie power (p) parameter is now in `print()` and `tidy()` output.

* `thetaf` is now `tweedie_p` in `sdmTMB_sim()`.

# sdmTMB 0.0.19.9003

* Fix bug affecting prediction with `se_fit = TRUE` for breakpoint models.

# sdmTMB 0.0.19.9002

* Simulation from the parameter covariance matrix works if random effects
  are turned off. #57

# sdmTMB 0.0.19.9000

* Smoothers `s()` are now *penalized* smoothers: they determine the 
  degree of wiggliness (as in mgcv) and it is no longer necessary to
  choose an appropriate `k` value a priori. Models fit with previous
  versions of sdmTMB with  `s(x, k = ...)` will not match models
  specified the same way in version >= 0.0.19 since the basis functions
  are now penalized. All the various `mgcv::s()` options should be supported
  but `t2()` (and `ti()` and `te()`) are not supported.

# sdmTMB 0.0.18.9001

* Add ELPD (expected log predictive density) to `sdmTMB_cv()`
  <https://arxiv.org/abs/1507.04544>
  
* Fix bug evaluating `...` when `sdmTMB_cv()` was called within a function. #54

# sdmTMB 0.0.18.9000

* Fix minor error in PC Matern prior

# sdmTMB 0.0.17.9000

* Add random walk option: `fields = "RW"`.

* Depreciate `ar1_fields` argument. See new `fields` argument in `sdmTMB().

* Many packages moved from 'Imports' to 'Suggests'

# sdmTMB 0.0.16.9000

* Lower default `nlminb()` `eval.max` and `iter.max` to 1000 and 2000.

* Added `profile` option in `sdmTMBcontrol()`. This can dramatically
  improve model fitting speed with many fixed effects. Note the
  result is likely to be slightly different with `TRUE` vs. `FALSE`.

* Added simulation from the MVN precision matrix to `predict.sdmTMB()`. 
  See the `sims` argument.
  
* Added `gather_sims()` and `spread_sims()` to extract parameter
  simulations from the joint precision matrix in a format that
  matches the tidybayes package.

* Added `get_index_sims()` for a population index calculated from
  the MVN simulation draws.
  
* Added `extract_mcmc()` to extract MCMC samples if the model is
  passed to tmbstan.
  
* Added the ability to predict from a model fitted with tmbstan.
  See the `tmbstan_model` argument in `predict.sdmTMB()`.

* Allowed for separate random field Matern range parameters for 
  spatial and spatiotemporal fields. E.g. `sdmTMB(shared_range = FALSE)`

* Bounded the AR1 rho parameter between -0.999 and 0.999 to improve 
  convergence; was -1 to 1. Please post an issue if this creates
  problems for your model.

* Added `map`, `start`, `lower`, and `upper` options to control
  model fitting. See `sdmTMBcontrol()`.

* Added priors for all parameters. See `?sdmTMB::priors` and the
  `priors` argument in `sdmTMB()`. PC priors are available for
  the random fields. See `?pc_matern` and the details there.
  
* Moved many less-common arguments from `sdmTMB()` to `sdmTMBcontrol()`.

* Fix bug in `sdmTMB_cv()` where fitting and testing data splits
  were reversed. I.e., the small chunk was fit; the big chunk was tested.

# sdmTMB 0.0.15.9000

* Added experimental penalized complexity (PC) prior as used in INLA.
  See arguments `matern_prior_O` and `matern_prior_E`.

* Added back `normalize` argument to `sdmTMB()` and default to `FALSE`.
  Setting to `TRUE` can dramatically speed up some model fits
  (~4 times for some test models).

# sdmTMB 0.0.14.9003

* Added vignette on making pretty maps of the output

# sdmTMB 0.0.14.9001

* Added some protections for possible user errors:
  * AR1 with a spatial-only model
  * Missing factor levels in time
  * Coordinate systems that are too big

# sdmTMB 0.0.14.9000

* Add `re_form_iid` to `predict.sdmTMB()`.

* Add `map_rf` option to `sdmTMB()`. This lets you map (fix at
  their starting values of zero) all random fields to produce a
  classic GLM/GLMM.

# sdmTMB 0.0.13.9000

* Add IID random intercepts interface. E.g. `... + (1 | g)` #34

# sdmTMB 0.0.12.9000

* Add `epsilon_predictor` argument in `sdmTMB()` to allow a model of the
  spatiotemporal variance through time.

# sdmTMB 0.0.11.9000

* Add `penalties` argument to allow for regularization.

# sdmTMB 0.0.10.9001

* Fix Student-t degrees of freedom in the randomized quantile residuals

# sdmTMB 0.0.10.9000

* Fixed parameter initialization for inverse links #35

* Switched Gamma 'phi' parameter to representing shape instead of CV to
  match glm(), glmmTMB(), etc.

# sdmTMB 0.0.9.9000

* Switched the density/abundance index calculation to use the link function as
  opposed to a hardcoded log() so that the `get_generic()` function can be used
  to grab things like standardized average values of the response across a grid.
  What used to be `log_total` in the raw TMB output is now `link_total` but most
  users you shouldn't notice any difference.

# sdmTMB 0.0.8.9000

* Overhauled the simulation function. The function is now called `sdmTMB_sim()`
  and uses INLA functions instead of RandomFields functions for simulating
  the random fields.

* The simulation function can now accommodate all families and links and takes
  an INLA mesh as input.

# sdmTMB 0.0.7.9001

* Allow specifying degrees of freedom in the Student-t family #29

# sdmTMB 0.0.7.9000

* Added a `tidy()` method (from broom and broom.mixed) to return a data frame
  of parameter estimates. The function can extract the fixed effects or the
  random effect parameters (variances, AR1 correlation, spatial range).

* Added an argument `extra_time` to `sdmTMB()`. This introduces additional time
  slices that you can then predict on if you want to interpolate or forecast.
  Internally, it uses Eric Ward's 'weights hack'. This is also useful if you
  have data unevenly spaced in time and you want the gaps evenly spaced for a
  random walk or AR1 process (add any missing years to `extra_time`).

* `make_spde()` is now replaced with `make_mesh()` and `make_spde()` has been
  soft deprecated. `make_mesh()` carries through the x and y column names to
  the predict function and is more in line with the tidyverse style of taking a
  data frame first.

* `make_mesh()` can accept `cutoff` as an argument (as in INLA), which is
  likely a better default way to specify the mesh since it scales across
  regions better and is line with the literature on INLA.

* `make_mesh()` can use a binary search algorithm to find a cutoff that best
  matches a desired number of knots (thanks to Kelli Johnson for the idea).

* Barrier meshes are now possible. See `add_barrier_mesh()` for an example.

* There is a pkgdown website now that gets auto generated with GitHub actions.

* There is the start of a model description vignette.
  It is very much a work in progress.

# sdmTMB 0.0.6.9009

* Fixed bug in dlnorm

# sdmTMB 0.0.6.9005

* Fixed bug in predictions with standard errors where one(?)
  parameter (a breakpoint parameter) would be passed in at its initial
  instead of MLE value.

# sdmTMB 0.0.6.9004

* Fixed bug with predictions on new data in models with break points

* Overhauled cross validation function. The function now:
    * uses Eric's weights hack so it can also be used for forecasting
    * initializes subsequent folds at the MLE of the first fold for
      considerable speed increases
    * works in parallel if a future plan initialized; see examples

* Added threshold parameters to the print method

* Added forecasting example with the weights hack

* Fixed bug in linear break point models

#  sdmTMB 0.0.6.9002

* Fixed GAM predictions with all 0s in new data.

* Add linear and logistic threshold models. #17

# sdmTMB 0.0.5.9000

* Added parsing of mgcv formulas for splines. #16

* Added ability to predict with standard errors at the population level. This
  helps with making marginal-effect plots. #15

* Added optimization options to aid convergence. Also added
  `run_extra_optimization()` to run these on already fit models. Default is for
  no extra optimization.

* Added binomial likelihood to cross validation. Git hash `ee3f3ba`.

* Started keeping track of news in `NEWS.md`.

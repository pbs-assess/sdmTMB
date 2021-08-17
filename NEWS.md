# sdmTMB

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

# sdmTMB 0.0.15.9000

* Added experimental penalized complexity (PC) prior as used in INLA.
  See arguments `matern_prior_O` and `matern_prior_E`.

* Added back `normalize` argument to `sdmTMB()` and default to `FALSE`.
  Setting to `TRUE` can dramatically speed up some model fits
  (~4 times for some test models).
>>>>>>> priors-experimental

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
  soft depreciated. `make_mesh()` carries through the x and y column names to
  the predict function and is more in line with the tidyverse style of taking a
  data frame first.

* `make_mesh()` can accept `cutoff` as an argument (as in INLA), which is
  likely a better default way to specify the mesh since it scales across
  regions better and is line with the literature on INLA.

* `make_mesh()` can use a binary search algorithm to find a cutoff that best
  matches a desired number of knots (thanks to Kelli Johnson for the idea).

* Barrier meshes are now possible. See `add_barrier_mesh()` for an example.

* There is a pkgdown website now that gets auto generated with GitHub actions:
  <https://pbs-assess.github.io/sdmTMB/index.html>

* There is the start of a model description vignette:
  <https://github.com/pbs-assess/sdmTMB/blob/devel/vignettes/model-description.Rmd>
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

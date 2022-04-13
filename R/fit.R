#' @useDynLib sdmTMB, .registration = TRUE
NULL

#' Fit a spatial or spatiotemporal GLMM with TMB
#'
#' Fit a spatial or spatiotemporal Gaussian random field GLMM with TMB using
#' the SPDE approach.
#' This can be useful for (dynamic) species distribution models and relative
#' abundance index standardization among many other uses.
#'
#' @param formula Model formula. See the Details section below for how to specify
#'   offsets and threshold parameters. For index standardization, you may wish
#'   to include `0 + as.factor(year)` (or whatever the time column is called)
#'   in the formula. IID random intercepts are possible using \pkg{lme4}
#'   syntax, e.g., `+ (1 | g)` where `g` is a column with factor levels.
#'   Penalized splines are possible via \pkg{mgcv} with `s()`. See examples
#'   and details below.
#' @param data A data frame.
#' @param mesh An object from [make_mesh()].
#' @param time An optional time column name (as character). Can be left as
#'   `NULL` for a model with only spatial random fields unless you wish to use
#'   one of the index or center of gravity functions over time.
#' @param family The family and link. Supports [gaussian()], [Gamma()],
#'   [binomial()], [poisson()], \code{\link[sdmTMB:families]{Beta()}},
#'   \code{\link[sdmTMB:families]{nbinom2()}},
#'   \code{\link[sdmTMB:families]{truncated_nbinom2()}},
#'   \code{\link[sdmTMB:families]{nbinom1()}},
#'   \code{\link[sdmTMB:families]{truncated_nbinom1()}},
#'   \code{\link[sdmTMB:families]{censored_poisson()}},
#'   \code{\link[sdmTMB:families]{student()}}, and
#'   \code{\link[sdmTMB:families]{tweedie()}}]. For binomial family options,
#'   see the 'Binomial families' in the Details section below.
#' @param spatial Estimate spatial random fields? Options are
#'   `'on'` / `'off'` or `TRUE` / `FALSE`.
#' @param spatiotemporal Estimate the spatiotemporal random fields as `'IID'`
#'   (independent and identically distributed; default), stationary `'AR1'`
#'   (first-order autoregressive), as a random walk (`'RW'`), or as fixed at 0
#'   `'off'`. Will be set to `'off'` if `time = NULL`. Note that the
#'   spatiotemporal standard deviation represents the marginal steady-state
#'   standard deviation of the process in the case of the AR1. I.e., it is
#'   scaled according to the correlation. See the [TMB
#'   documentation](https://kaskr.github.io/adcomp/classAR1__t.html). If the AR1
#'   correlation coefficient (rho) is estimated close to 1, say > 0.99, then you
#'   may wish to switch to the random walk `"RW"`. Capitalization is ignored. `TRUE`
#'   gets converted to `'iid'` and `FALSE` gets converted to `off`.
#' @param share_range Logical: estimate a shared spatial and spatiotemporal
#'   range parameter (`TRUE`) or independent range parameters (`FALSE`).
#' @param time_varying An optional one-sided formula describing covariates that
#'   should be modelled as a random walk through time. Be careful not to include
#'   covariates (including the intercept) in both the main and time-varying
#'   formula. I.e., at least one should have `~ 0` or `~ -1`.
#' @param spatial_varying An optional one-sided formula **with a single
#'   predictor** of a coefficient that should varying in space as a random
#'   field. Note that you may want to include a fixed effect for the same
#'   variable to improve interpretability. If the (scaled) time column, will
#'   represent a local-time-trend model. See \doi{10.1111/ecog.05176} and the
#'   [spatial trends
#'   vignette](https://pbs-assess.github.io/sdmTMB/articles/spatial-trend-models.html).
#'    Note this predictor should be centered to have mean zero and have a
#'   standard deviation of approximately 1 (scale by the SD).
#' @param weights A numeric vector representing optional likelihood weights for
#'   the conditional model. Implemented as in \pkg{glmmTMB}: weights do not have
#'   to sum to one and are not internally modified. Can also be used for trials
#'   with the binomial family; the weights argument needs to be a vector and not
#'   a name of the variable in the data frame. See the Details section below.
#' @param offset A numeric vector representing the model offset. In delta/hurdle
#'   models, this applies only to the positive component. Usually a log
#'   transformed variable. Not used in any prediction.
#' @param extra_time Optional extra time slices (e.g., years) to include for
#'   interpolation or forecasting with the predict function. See the
#'   Details section below.
#' @param reml Logical: use REML (restricted maximum likelihood) estimation
#'   rather than maximum likelihood? Internally, this adds the fixed effects
#'   to the list of random effects to integrate over.
#' @param silent Silent or include optimization details?
#' @param anisotropy Logical: allow for anisotropy? See [plot_anisotropy()].
#' @param control Optimization control options via [sdmTMBcontrol()].
#' @param priors Optional penalties/priors via [sdmTMBpriors()].
#' @param previous_fit A previously fitted sdmTMB model to initialize the
#'   optimization with. Can greatly speed up fitting. Note that the data and
#'   model must be set up *exactly* the same way. However, the `weights` argument
#'   can change, which can be useful for cross-validation.
#' @param do_fit Fit the model (`TRUE`) or return the processed data without
#'   fitting (`FALSE`)?
#' @param experimental A named list for esoteric or in-development options.
#'    Here be dragons.
#'   (Experimental) A column name (as character) of a
#'   predictor of a linear trend (in log space) of the spatiotemporal standard
#'   deviation. By default, this is `NULL` and fits a model with a constant
#'   spatiotemporal variance. However, this argument can also be a character
#'   name in the original data frame (a covariate that ideally has been
#'   standardized to have mean 0 and standard deviation = 1). Because the
#'   spatiotemporal field varies by time step, the standardization should be
#'   done by time. If the name of a predictor is included, a log-linear model is
#'   fit where the predictor is used to model effects on the standard deviation,
#'   e.g. `log(sd(i)) = B0 + B1 * epsilon_predictor(i)`. The 'epsilon_model' argument may also
#' be specified. This is the name of the model to use to modeling time-varying epsilon. This
#' can be one of the following: "trend" (default, fits a linear model without random effects),
#' "re" (fits a model with random effects in epsilon_st, but no trend), and "trend-re" (a model
#' that includes both the trend and random effects)
#' @param fields **Depreciated.** Replaced by `spatiotemporal`.
#' @param include_spatial **Depreciated.** Replaced by `spatial`.
#' @param spatial_only **Depreciated.** Replaced by `spatiotemporal = "off"`.
#' @param spde **Depreciated.** Replaced by `mesh`.
#' @param ... Not currently used.
#' @importFrom methods as is
#' @importFrom mgcv s t2
#' @importFrom stats gaussian model.frame model.matrix as.formula
#'   model.response terms model.offset
#' @importFrom lifecycle deprecated is_present deprecate_warn
#'
#' @return
#' An object (list) of class `sdmTMB`. Useful elements include:
#'
#' * `sd_report`: output from [TMB::sdreport()]
#' * `gradients`: log likelihood gradients with respect to each fixed effect
#' * `model`: output from [stats::nlminb()]
#' * `data`: the fitted data
#' * `mesh`: the object that was supplied to the `mesh` argmument
#' * `family`: the family object, which includes the inverse link function
#' * `tmb_params`: The parameters list passed to [TMB::MakeADFun()]
#' * `tmb_map`: The 'map' list passed to [TMB::MakeADFun()]
#' * `tmb_data`: The data list passed to [TMB::MakeADFun()]
#' * `tmb_obj`: The TMB object created by [TMB::MakeADFun()]
#'
#' @details
#'
#' **Model description**
#'
#' For now, see the
#' [model description](https://pbs-assess.github.io/sdmTMB/articles/model-description.html)
#' vignette for a start. There are also descriptions of particular models in
#' Anderson et al. (2019) and Barnett et al. (2020) (see reference list below).
#'
#' **Binomial families**
#'
#' Following the structure of [stats::glm()] and \pkg{glmmTMB}, a binomial
#' family can be specified in one of 4 ways: (1) the response may be a factor
#' (and the model classifies the first level versus all others), (2) the
#' response may be binomial (0/1), (3) the response can be a matrix of form
#' `cbind(success, failure)`, and (4) the response may be the observed
#' proportions, and the 'weights' argument is used to specify the Binomial size
#' (N) parameter (`prob ~ ..., weights = N`).
#'
#' **Smooth terms**
#'
#' Smooth terms can be included following GAMs (generalized additive models)
#' using `+ s(x)`, which implements a smooth from [mgcv::s()]. \pkg{sdmTMB} uses
#' penalized smooths, constructed via [mgcv::smooth2random()]. This is a similar
#' approach implemented in \pkg{gamm4} and \pkg{brms}, among other packages.
#' Within these smooths, the same syntax commonly used in [mgcv::s()] can be
#' applied, e.g. 2-dimensional smooths may be constructed with `+ s(x, y)`;
#' smooths can be specific to various factor levels, `+ s(x, by = group)`; the
#' basis function dimensions may be specified, e.g. `+ s(x, k = 4)`; and various
#' types of splines may be constructed such as cyclic splines to model
#' seasonality, `+ s(month, bs = "cc", k = 12)`. Prior to version 0.0.18,
#' \pkg{sdmTMB} implemented unpenalized splines.
#'
#' **Threshold models**
#'
#' A linear break-point relationship for a covariate can be included via
#' `+ breakpt(variable)` in the formula, where `variable` is a single covariate
#' corresponding to a column in `data`. In this case, the relationship is linear
#' up to a point and then constant (hockey-stick shaped).
#'
#' Similarly, a logistic-function threshold model can be included via
#' `+ logistic(variable)`. This option models the relationship as a logistic
#' function of the 50% and 95% values. This is similar to length- or size-based
#' selectivity in fisheries, and is parameterized by the points at which f(x) =
#' 0.5 or 0.95. See the vignette.
#'
#' Note that only a single threshold covariate can be included.
#'
#' See the
#' [threshold vignette](https://pbs-assess.github.io/sdmTMB/articles/threshold-models.html).
#'
#' **Extra time: forecasting or interpolating**
#'
#' Extra time slices (e.g., years) can be included for interpolation or
#' forecasting with the predict function via the `extra_time` argument. The
#' predict function requires all time slices to be defined when fitting the
#' model to ensure the various time indices are set up correctly. Be careful if
#' including extra time slices that the model remains identifiable. For example,
#' including `+ as.factor(year)` in `formula` will render a model with no data
#' to inform the expected value in a missing year. [sdmTMB()] makes no attempt
#' to determine if the model makes sense for forecasting or interpolation. The
#' options `time_varying`, `spatiotemporal = "RW"`, and `spatiotemporal = "AR1"`
#' provide mechanisms to predict over missing time slices with process error.
#'
#' `extra_time` can also be used to fill in missing time steps for the purposes
#' of a random walk or AR1 spatiotemporal field if their inclusion makes the gaps
#' between time steps even.
#'
#' **Index standardization**
#'
#' For index standardization, you may wish to include `0 + as.factor(year)`
#' (or whatever the time column is called) in the formula. See a basic
#' example of index standardization in the relevant
#' [package vignette](https://pbs-assess.github.io/sdmTMB/articles/model-description.html).
#' You will need to specify the `time` argument. See [get_index()] and/or
#' [get_index_sims()].
#'
#' **Regularization and priors**
#'
#' You can achieve regularization via penalties (priors) on the fixed effect
#' parameters. See [sdmTMBpriors()]. These should not include `offset` terms.
#' You can fit the model once without penalties and inspect the element
#' `head(your_model$tmb_data$X_ij)` if you want to see how the formula is
#' translated to the fixed effect model matrix.
#'
#' @references
#'
#' Main reference introducing the package to cite when using sdmTMB:
#'
#' Anderson, S.C., E.J. Ward, P.A. English, L.A.K. Barnett. 2022. sdmTMB: an R
#' package for fast, flexible, and user-friendly generalized linear mixed effects
#' models with spatial and spatiotemporal random fields.
#' bioRxiv 2022.03.24.485545; \doi{10.1101/2022.03.24.485545}.
#'
#' Reference for local trends:
#'
#' Barnett, L.A.K., E.J. Ward, S.C. Anderson. Improving estimates of species
#' distribution change by incorporating local trends. Ecography. 44(3):427-439.
#' \doi{10.1111/ecog.05176}.
#'
#' Further explanation of the model and application to calculating climate
#' velocities:
#'
#' English, P., E.J. Ward, C.N. Rooper, R.E. Forrest, L.A. Rogers, K.L. Hunter,
#' A.M. Edwards, B.M. Connors, S.C. Anderson. 2021. Contrasting climate velocity
#' impacts in warm and cool locations show that effects of marine warming are
#' worse in already warmer temperate waters. In press at Fish and Fisheries.
#' \doi{10.1111/faf.12613}.
#'
#' Code for implementing the barrier-SPDE written by Olav Nikolai Breivik and
#' Hans Skaug.
#'
#' A number of sections of the original TMB model code were adapted from the
#' VAST R package:
#'
#' Thorson, J.T., 2019. Guidance for decisions using the Vector Autoregressive
#' Spatio-Temporal (VAST) package in stock, ecosystem, habitat and climate
#' assessments. Fish. Res. 210:143–161.
#' \doi{10.1016/j.fishres.2018.10.013}.
#'
#' Code for the `family` R-to-TMB implementation, selected parameterizations of
#' the observation distributions, general package structure inspiration, and the
#' idea behind the TMB prediction approach were adapted from the glmmTMB R
#' package:
#'
#' Mollie E. Brooks, Kasper Kristensen, Koen J. van Benthem, Arni Magnusson,
#' Casper W. Berg, Anders Nielsen, Hans J. Skaug, Martin Maechler and Benjamin
#' M. Bolker (2017). glmmTMB Balances Speed and Flexibility Among Packages for
#' Zero-inflated Generalized Linear Mixed Modeling. The R Journal, 9(2):378-400.
#' \doi{10.32614/rj-2017-066}.
#'
#' @export
#'
#' @examples
#' if (inla_installed()) {
#'
#' # Build a fairly coarse mesh for example speed:
#' mesh <- make_mesh(pcod_2011, c("X", "Y"), cutoff = 20)
#' # `cutoff = 10` or `15` might make more sense in applied situations
#' # `cutoff` is the minimum distance between mesh vertices in units of the
#' # x and y coordinates.
#' # Or build any mesh in INLA and pass it to the `mesh` argument.
#'
#' # Quick mesh plot:
#' plot(mesh)
#' # Or:
#' # ggplot2::ggplot() + inlabru::gg(mesh$mesh)
#'
#' # Fit a Tweedie spatial random field GLMM with a smoother for depth:
#' fit <- sdmTMB(
#'   density ~ s(depth),
#'   data = pcod_2011, mesh = mesh,
#'   family = tweedie(link = "log")
#' )
#' fit
#'
#' # Extract coefficients:
#' tidy(fit, conf.int = TRUE)
#' tidy(fit, effects = "ran_par", conf.int = TRUE)
#'
#' # Check maximum gradient is small:
#' max(fit$gradients)
#'
#' # Predict; see ?predict.sdmTMB
#' p <- predict(fit)
#' p <- predict(fit, newdata = subset(qcs_grid, year == 2017))
#'
#' # Add spatiotemporal random fields:
#' fit <- sdmTMB(
#'   density ~ 1 + s(depth) + as.factor(year),
#'   time = "year", #<
#'   data = pcod_2011, mesh = mesh,
#'   family = tweedie(link = "log")
#' )
#' fit
#'
#' # Make the fields AR1:
#' fit <- sdmTMB(
#'   density ~ s(depth),
#'   time = "year",
#'   spatial = "off",
#'   spatiotemporal = "ar1", #<
#'   data = pcod_2011, mesh = mesh,
#'   family = tweedie(link = "log")
#' )
#' fit
#'
#' # Make the fields a random walk:
#' fit <- sdmTMB(
#'   density ~ s(depth),
#'   time = "year",
#'   spatial = "off",
#'   spatiotemporal = "rw", #<
#'   data = pcod_2011, mesh = mesh,
#'   family = tweedie(link = "log")
#' )
#' fit
#'
#' # Depth smoothers by year:
#' fit <- sdmTMB(
#'   density ~ s(depth, by = as.factor(year)), #<
#'   time = "year",
#'   spatial = "off",
#'   spatiotemporal = "rw",
#'   data = pcod_2011, mesh = mesh,
#'   family = tweedie(link = "log")
#' )
#' fit
#'
#' # 2D depth-year smoother:
#' fit <- sdmTMB(
#'   density ~ s(depth, year), #<
#'   spatial = "off",
#'   data = pcod_2011, mesh = mesh,
#'   family = tweedie(link = "log")
#' )
#' fit
#'
#' # Turn off spatial random fields:
#' fit <- sdmTMB(
#'   present ~ poly(log(depth)),
#'   spatial = "off", #<
#'   data = pcod_2011, mesh = mesh,
#'   family = binomial()
#' )
#' fit
#'
#' # Which, matches glm():
#' fit_glm <- glm(
#'   present ~ poly(log(depth)),
#'   data = pcod_2011,
#'   family = binomial()
#' )
#' summary(fit_glm)
#' AIC(fit, fit_glm)
#'
#' # Delta/hurdle binomial-Gamma model:
#' fit_dg <- sdmTMB(
#'   density ~ poly(log(depth)),
#'   data = pcod_2011, mesh = mesh,
#'   spatial = "off",
#'   family = delta_gamma() #<
#' )
#' fit_dg
#'
#' # Delta/hurdle Poisson-link (cloglog link):
#' fit_pg <- sdmTMB(
#'   density ~ s(depth),
#'   data = pcod_2011, mesh = mesh,
#'   spatial = "off",
#'   family = delta_poisson_link_gamma() #<
#' )
#' fit_pg
#'
#' # Delta/hurdle truncated NB2:
#' pcod_2011$count <- round(pcod_2011$density)
#' fit_nb2 <- sdmTMB(
#'   count ~ s(depth),
#'   data = pcod_2011, mesh = mesh,
#'   spatial = "off",
#'   family = delta_truncated_nbinom2() #<
#' )
#' fit_nb2
#'
#' # Regular NB2:
#' fit_nb2 <- sdmTMB(
#'   count ~ s(depth),
#'   data = pcod_2011, mesh = mesh,
#'   spatial = "off",
#'   family = nbinom2() #<
#' )
#' fit_nb2
#'
#' # IID random intercepts by year:
#' pcod_2011$fyear <- as.factor(pcod_2011$year)
#' fit <- sdmTMB(
#'   density ~ s(depth) + (1 | fyear), #<
#'   data = pcod_2011, mesh = mesh,
#'   family = tweedie(link = "log")
#' )
#' fit
#'
#' # Spatially varying coefficient of year:
#' pcod_2011$year_scaled <- as.numeric(scale(pcod_2011$year))
#' fit <- sdmTMB(
#'   density ~ year_scaled,
#'   spatial_varying = ~ 0 + year_scaled, #<
#'   data = pcod_2011, mesh = mesh, family = tweedie(), time = "year"
#' )
#' fit
#'
#' # Time-varying effects of depth and depth squared:
#' fit <- sdmTMB(
#'   density ~ 0 + as.factor(year),
#'   time_varying = ~ 0 + depth_scaled + depth_scaled2, #<
#'   data = pcod_2011, time = "year", mesh = mesh,
#'   family = tweedie()
#' )
#' print(fit)
#' # Extract values:
#' est <- as.list(fit$sd_report, "Estimate")
#' se <- as.list(fit$sd_report, "Std. Error")
#' est$b_rw_t[, , 1]
#' se$b_rw_t[, , 1]
#'
#' # Linear break-point effect of depth:
#' fit <- sdmTMB(
#'   density ~ breakpt(depth_scaled), #<
#'   data = pcod_2011,
#'   # spatial = "off",
#'   mesh = mesh,
#'   family = tweedie()
#' )
#' fit
#' }

sdmTMB <- function(
  formula,
  data,
  mesh,
  time = NULL,
  family = gaussian(link = "identity"),
  spatial = c("on", "off"),
  spatiotemporal = c("IID", "AR1", "RW", "off"),
  share_range = TRUE,
  time_varying = NULL,
  spatial_varying = NULL,
  weights = NULL,
  offset = NULL,
  extra_time = NULL,
  reml = FALSE,
  silent = TRUE,
  anisotropy = FALSE,
  control = sdmTMBcontrol(),
  priors = sdmTMBpriors(),
  previous_fit = NULL,
  experimental = NULL,
  do_fit = TRUE,
  spatial_only = deprecated(),
  fields = deprecated(),
  include_spatial = deprecated(),
  spde = deprecated(),
  ...
  ) {


  delta <- isTRUE(family$delta)
  n_m <- if (isTRUE(delta)) 2L else 1L

  check_spatiotemporal_arg <- function(x, .which = 1) {
    sp_len <- length(x)
    .spatiotemporal <- tolower(as.character(x[[.which]]))
    assert_that(.spatiotemporal %in% c("iid", "ar1", "rw", "off", "true", "false"),
      msg = "`spatiotemporal` argument value not valid")
    if (.spatiotemporal == "true") .spatiotemporal <- "iid"
    if (.spatiotemporal == "false") .spatiotemporal <- "off"
    .spatiotemporal <- match.arg(tolower(.spatiotemporal), choices = c("iid", "ar1", "rw", "off"))
    if (is.null(time) && .spatiotemporal != "off" && sp_len >= 1) {
      stop("`spatiotemporal` is set but the `time` argument is missing.", call. = FALSE)
    }
    .spatiotemporal
  }

  if (!missing(spatiotemporal)) {
    if (delta && !is.list(spatiotemporal)) {
      spatiotemporal <- rep(spatiotemporal[[1]], 2L)
    }
    spatiotemporal <- vapply(seq_along(spatiotemporal),
      function(i) check_spatiotemporal_arg(spatiotemporal, i), FUN.VALUE = character(1L))
  } else {
    if (is.null(time))
      spatiotemporal <- rep("off", n_m)
    else
      spatiotemporal <- rep("iid", n_m)
  }

  # if (is_present(spatial_only)) {
  #   deprecate_warn("0.0.20", "sdmTMB(spatial_only)", "sdmTMB(spatiotemporal)", details = "`spatiotemporal = 'off'` (or `time = NULL`) is equivalent to `spatial_only = TRUE`.")
  #} else {

  if (is.null(time)) {
    spatial_only <- rep(TRUE, n_m)
  } else {
    spatial_only <- ifelse(spatiotemporal == "off", TRUE, FALSE)
  }
  # }
  # if (is_present(fields)) {
  #   deprecate_warn("0.0.20", "sdmTMB(fields)", "sdmTMB(spatiotemporal)")
  #   spatiotemporal <- tolower(fields)
  # }
  # if (is_present(include_spatial)) {
  #   deprecate_warn("0.0.20", "sdmTMB(include_spatial)", "sdmTMB(spatial)")
  #} else {
    if (!is.logical(spatial[[1]])) spatial <- match.arg(tolower(spatial), choices = c("on", "off"))
    if (identical(spatial, "on") || isTRUE(spatial)) {
      include_spatial <- TRUE
    } else {
      include_spatial <- FALSE
    }
  # }
  if (!include_spatial && all(spatiotemporal == "off") || !include_spatial && all(spatial_only)) {
    # message("Both spatial and spatiotemporal fields are set to 'off'.")
    control$map_rf <- TRUE
    if (missing(mesh)) {
      data$sdmTMB_X_ <- data$sdmTMB_Y_ <- stats::runif(nrow(data))
      mesh <- make_mesh(data, c("sdmTMB_X_", "sdmTMB_Y_"), cutoff = 1)
    }
  }

  share_range <- unlist(share_range)
  if (length(share_range) == 1L) share_range <- rep(share_range, n_m)
  share_range[spatiotemporal == "off"] <- TRUE

  if (is_present(spde)) {
    deprecate_warn("0.0.20", "sdmTMB(spde)", "sdmTMB(mesh)")
  } else {
    spde <- mesh
  }
  epsilon_model <- NULL
  epsilon_predictor <- NULL
  if (!is.null(experimental)) {
    if ("epsilon_predictor" %in% names(experimental)) {
      epsilon_predictor <- experimental$epsilon_predictor
    } else {
      epsilon_predictor <- NULL
    }

    if ("epsilon_model" %in% names(experimental)) {
      epsilon_model <- experimental$epsilon_model
    } else {
      epsilon_model <- NULL
    }

    if ("lwr" %in% names(experimental) && "upr" %in% names(experimental)) {
      lwr <- experimental$lwr
      upr <- experimental$upr
    } else {
      #epsilon_predictor <- NULL
      lwr <- 0
      upr <- Inf
    }
  } else {
    #epsilon_predictor <- NULL
    lwr <- 0
    upr <- Inf
  }
  normalize <- control$normalize
  nlminb_loops <- control$nlminb_loops
  newton_loops <- control$newton_loops
  quadratic_roots <- control$quadratic_roots
  start <- control$start
  multiphase <- control$multiphase
  map_rf <- control$map_rf
  map <- control$map
  lower <- control$lower
  upper <- control$upper
  get_joint_precision <- control$get_joint_precision
  dots <- list(...)

  if ("ar1_fields" %in% names(dots)) {
    lifecycle::deprecate_stop("0.0.20", "sdmTMB(ar1_fields)", "sdmTMB(spatiotemporal)")
  }
  if ("penalties" %in% names(dots)) {
    stop("`penalties` are now specified via the `priors` argument.",
      "E.g., `priors = sdmTMBpriors(b = normal(c(0, 0), c(1, 1)))`",
      "for 2 fixed effects.", call. = FALSE)
  }
  dot_checks <- c("lower", "upper", "profile", "parallel",
    "nlminb_loops", "newton_steps", "mgcv", "quadratic_roots", "multiphase",
    "newton_loops", "start", "map", "map_rf", "get_joint_precision", "normalize")
  .control <- control
  for (i in dot_checks) .control[[i]] <- NULL # what's left should be for nlminb
  dot_old <- dot_checks %in% names(dots)
  if (any(dot_old)) {
    stop("The ", if (sum(dot_old) > 1) "arguments" else "argument", " `",
      paste(dot_checks[dot_old], collapse = "`, `"),
      "` ", if (sum(dot_old) > 1) "were" else "was",
      " found in the call to `sdmTMB()`.\n",
      if (sum(dot_old) > 1) "These are " else "This is ",
      "now passed through the `control` argument via the\n",
      "`sdmTMBcontrol()` list.", call. = FALSE)
  }

  ar1_fields <- spatiotemporal == "ar1"
  rw_fields <- spatiotemporal == "rw"
  assert_that(
    is.logical(reml), is.logical(anisotropy), is.logical(share_range), is.logical(silent),
    is.logical(silent),
    is.logical(multiphase),
    is.logical(map_rf), is.logical(normalize)
  )
  if (!is.null(spatial_varying)) assert_that(class(spatial_varying) == "formula")
  assert_that(is.list(priors))
  if (!is.null(time_varying)) assert_that(identical(class(time_varying), "formula"))
  assert_that(is.list(.control))
  if (!is.null(previous_fit)) assert_that(identical(class(previous_fit), "sdmTMB"))
  if (!is.null(time)) assert_that(is.character(time))
  assert_that(identical(class(spde), "sdmTMBmesh"))
  assert_that(identical(class(formula), "formula"))
  assert_that("data.frame" %in% class(data))
  if (!is.null(map) && length(map) != length(start)) {
    warning("`length(map) != length(start)`. You likely want to specify ",
      "`start` values if you are setting the `map` argument.", call. = FALSE)
  }
  if (family$family[1] == "censored_poisson") {
    assert_that("lwr" %in% names(experimental) && "upr" %in% names(experimental),
      msg = "`lwr` and `upr` must be specified in `experimental` as elements of a named list to use the censored Poisson likelihood.")
    assert_that(length(lwr) == nrow(data) && length(upr) == nrow(data))
    assert_that(length(lwr) == length(upr))
    assert_that(mean(upr-lwr, na.rm = TRUE)>=0)
  }

  if (!is.null(time)) {
    assert_that(time %in% names(data),
      msg = "Specified `time` column is missing from `data`.")
  }
  if (is.null(time)) {
    time <- "_sdmTMB_time"
    data[[time]] <- 0L
  } else {
    if (sum(is.na(data[[time]])) > 1)
      stop("There is at least one NA value in the time column. ",
        "Please remove it.", call. = FALSE)
  }
  if (is.factor(data[[time]])) {
    if (length(levels(data[[time]])) > length(unique(data[[time]]))) {
      stop("The time column is a factor and there are extra factor levels.",
        "Please remove these or turn your time column into an integer.", call. = FALSE)
    }
  }
  assert_that(identical(nrow(spde$loc_xy), nrow(data)),
    msg = "Number of x-y coordinates in `mesh` does not match `nrow(data)`.")

  n_orig <- suppressWarnings(TMB::openmp(NULL))
  if (n_orig > 0 && .Platform$OS.type == "unix") { # openMP is supported
    TMB::openmp(n = control$parallel)
    on.exit({TMB::openmp(n = n_orig)})
  }

  thresh <- check_and_parse_thresh_params(formula, data)
  original_formula <- formula
  formula <- thresh$formula

  if (!is.null(extra_time)) { # for forecasting or interpolating
    if (!"xy_cols" %in% names(spde)) {
      stop("Please use make_mesh() instead of the deprecated make_mesh() to use `extra_time`.",
        call. = FALSE)
    }
    data <- expand_time(df = data, time_slices = extra_time, time_column = time)
    weights <- data$weight_sdmTMB
    spde$loc_xy <- as.matrix(data[,spde$xy_cols,drop=FALSE])
    spde$A_st <- INLA::inla.spde.make.A(spde$mesh, loc = spde$loc_xy)
    spde$sdm_spatial_id <- seq(1, nrow(data)) # FIXME
  }

  spatial_varying_formula <- spatial_varying # save it
  if (!is.null(spatial_varying)) {
    z_i <- model.matrix(spatial_varying, data)
    .int <- grep("(Intercept)", colnames(z_i))
    if (sum(.int) > 0) z_i <- z_i[,-.int,drop=FALSE]
    spatial_varying <- colnames(z_i)
  } else {
    z_i <- matrix(0, nrow(data), 0L)
  }
  n_z <- ncol(z_i)

  contains_offset <- check_offset(formula)
  if (contains_offset) warning("Contains offset in formula. This is deprecated. Please use the `offset` argument.", call. = FALSE)

  # Parse random intercepts:
  split_formula <- glmmTMB::splitForm(formula)
  RE_names <- barnames(split_formula$reTrmFormulas)
  fct_check <- vapply(RE_names, function(x) check_valid_factor_levels(data[[x]], .name = x), TRUE)
  RE_indexes <- vapply(RE_names, function(x) as.integer(data[[x]]) - 1L, rep(1L, nrow(data)))
  nobs_RE <- unname(apply(RE_indexes, 2L, max)) + 1L
  if (length(nobs_RE) == 0L) nobs_RE <- 0L
  formula <- split_formula$fixedFormula
  ln_tau_G_index <- unlist(lapply(seq_along(nobs_RE), function(i) rep(i, each = nobs_RE[i]))) - 1L

  formula_no_sm <- remove_s_and_t2(formula)
  X_ij <- model.matrix(formula_no_sm, data)
  mf <- model.frame(formula_no_sm, data)
  mt <- attr(mf, "terms")
  # parse everything mgcv + smoothers:
  sm <- parse_smoothers(formula = formula, data = data)

  # deprecated:
  offset_pos <- grep("^offset$", colnames(X_ij))
  y_i <- model.response(mf, "numeric")
  if (family$family[1] %in% c("Gamma", "lognormal") && min(y_i) <= 0 && !delta) {
    stop("Gamma and lognormal must have response values > 0.", call. = FALSE)
  }

  # This is taken from approach in glmmTMB to match how they handle binomial
  # yobs could be a factor -> treat as binary following glm
  # yobs could be cbind(success, failure)
  # yobs could be binary
  # (yobs, weights) could be (proportions, size)
  # On the C++ side 'yobs' must be the number of successes.
  size <- rep(1, nrow(X_ij)) # for non-binomial case
  if (identical(family$family[1], "binomial") && !delta) {
    ## call this to catch the factor / matrix cases
    y_i <- model.response(mf, type = "any")
    if (is.factor(y_i)) {
      ## following glm, ‘success’ is interpreted as the factor not
      ## having the first level (and hence usually of having the
      ## second level).
      y_i <- pmin(as.numeric(y_i) - 1, 1)
      size <- rep(1, length(y_i))
    } else {
      if (is.matrix(y_i)) { # yobs=cbind(success, failure)
        size <- y_i[, 1] + y_i[, 2]
        yobs <- y_i[, 1] # successes
        y_i <- yobs
      } else {
        if (all(y_i %in% c(0, 1))) { # binary
          size <- rep(1, length(y_i))
        } else { # proportions
          y_i <- weights * y_i
          size <- weights
          weights <- rep(1, length(y_i))
        }
      }
    }
  }

  if (identical(family$link[1], "log") && min(y_i, na.rm = TRUE) < 0 && !delta) {
    stop("`link = 'log'` but the reponse data include values < 0.", call. = FALSE)
  }

  if (is.null(offset)) offset <- rep(0, length(y_i))

  if (!is.null(time_varying)) {
    X_rw_ik <- model.matrix(time_varying, data)
  } else {
    X_rw_ik <- matrix(0, nrow = nrow(data), ncol = 1)
  }

  n_s <- nrow(spde$mesh$loc)

  barrier <- "spde_barrier" %in% names(spde)
  if (barrier && anisotropy) {
    warning("Using a barrier mesh; therefore, anistropy will be disabled.", call. = FALSE)
    anisotropy <- FALSE
  }
  df <- if (family$family[1] == "student" && "df" %in% names(family)) family$df else 3

  est_epsilon_model <- 0L
  epsilon_covariate <- rep(0, length(unique(data[[time]])))
  if (!is.null(epsilon_predictor) & !is.null(epsilon_model)) {
    if(epsilon_model %in% c("trend","trend-re")) {
      # covariate vector dimensioned by number of time steps
      time_steps <- unique(data[[time]])
      for (i in seq_along(time_steps)) {
        epsilon_covariate[i] <- data[data[[time]] == time_steps[i],
                                     epsilon_predictor, drop = TRUE][[1]]
      }
      est_epsilon_model <- 1L
    }
  }
  # flags for turning off the trend and random effects
  est_epsilon_slope <- 0
  if(!is.null(epsilon_model)) {
    if(epsilon_model %in% c("trend","trend-re")) {
      est_epsilon_slope <- 1L
      est_epsilon_model <- 1L
    }
  }
  est_epsilon_re <- 0
  if(!is.null(epsilon_model)) {
    if(epsilon_model[1] %in% c("re","trend-re")) {
      est_epsilon_re <- 1L
      est_epsilon_model <- 1L
    }
  }


  priors_b <- priors$b
  .priors <- priors
  .priors$b <- NULL # removes this in the list, so not passed in as data
  if (nrow(priors_b) == 1L && ncol(X_ij) > 1L) {
    if (!is.na(priors_b[[1]])) {
      message("Expanding `b` priors to match model matrix.")
    }
    # creates matrix that is 2 columns of NAs, rows = number of unique bs
    # Instead of passing in a 2-column matrix of NAs, pass in a matrix that
    # has means in first col and the remainder is Var-cov matrix
    priors_b <- mvnormal(rep(NA, ncol(X_ij)))
  }
  # ncol(X_ij) may occur if time varying model, no intercept
  if (ncol(X_ij) > 0 & !identical(nrow(priors_b), ncol(X_ij)))
    stop("The number of 'b' priors does not match the model matrix.", call. = FALSE)
  if (ncol(priors_b) == 2 && attributes(priors_b)$dist == "normal") {
    # normal priors passed in by user; change to MVN diagonal matrix
    # as.numeric() here prevents diag(NA) taking on factor
    if (length(priors_b[, 2]) == 1L) {
      if (is.na(priors_b[, 2])) priors_b[, 2] <- 1
    }
    priors_b <- mvnormal(location = priors_b[, 1], scale = diag(as.numeric(priors_b[, 2]), ncol = nrow(priors_b)))
  }
  # in some cases, priors_b will be a mix of NAs (no prior) and numeric values
  # easiest way to deal with this is to subset the Sigma matrix
  not_na <- which(!is.na(priors_b[,1]))
  if (length(not_na) == 0L) {
    priors_b_Sigma <- diag(2) # TMB can't handle 0-dim matrix
  } else {
    Sigma <- as.matrix(priors_b[, -1])
    priors_b_Sigma <- as.matrix(Sigma[not_na, not_na])
  }

  if (!"A_st" %in% names(spde)) stop("`mesh` was created with an old version of `make_mesh()`.", call. = FALSE)
  if (delta) y_i <- cbind(ifelse(y_i > 0, 1, 0), ifelse(y_i > 0, y_i, NA_real_))
  if (!delta) y_i <- matrix(y_i, ncol = 1L)

  tmb_data <- list(
    y_i        = y_i,
    n_t        = length(unique(data[[time]])),
    z_i        = z_i,
    offset_i   = offset,
    proj_offset_i = 0,
    A_st       = spde$A_st,
    sim_re     = if ("sim_re" %in% names(experimental)) as.integer(experimental$sim_re) else rep(0L, 6),
    A_spatial_index = spde$sdm_spatial_id - 1L,
    year_i     = make_year_i(data[[time]]),
    ar1_fields = ar1_fields,
    rw_fields =  rw_fields,
    X_ij       = X_ij,
    X_rw_ik    = X_rw_ik,
    Zs         = sm$Zs, # optional smoother basis function matrices
    Xs         = sm$Xs, # optional smoother linear effect matrix
    proj_Zs    = list(),
    proj_Xs    = matrix(nrow = 0L, ncol = 0L),
    b_smooth_start = sm$b_smooth_start,
    proj_lon   = 0,
    proj_lat   = 0,
    do_predict = 0L,
    calc_se    = 0L,
    pop_pred   = 0L,
    exclude_RE = rep(0L, ncol(RE_indexes)),
    weights_i  = if (!is.null(weights)) weights else rep(1, length(y_i)),
    area_i     = rep(1, length(y_i)),
    normalize_in_r = 0L, # not used first time
    flag = 1L, # part of TMB::normalize()
    calc_index_totals = 0L,
    calc_cog = 0L,
    random_walk = as.integer(!is.null(time_varying)),
    priors_b_n = length(not_na),
    priors_b_index = not_na - 1L,
    priors_b_mean = priors_b[not_na,1],
    priors_b_Sigma = priors_b_Sigma,
    priors = as.numeric(unlist(.priors)),
    share_range = as.integer(share_range),
    include_spatial = as.integer(include_spatial),
    proj_mesh  = Matrix::Matrix(0, 1, 1, doDiag = FALSE), # dummy
    proj_X_ij  = matrix(0, ncol = 1, nrow = 1), # dummy
    proj_X_rw_ik = matrix(0, ncol = 1, nrow = 1), # dummy
    proj_year  = 0, # dummy
    proj_spatial_index = 0, # dummy
    proj_z_i = matrix(0, nrow = 1, ncol = n_m), # dummy
    spde_aniso = make_anisotropy_spde(spde, anisotropy),
    spde       = spde$spde$param.inla[c("M0","M1","M2")],
    barrier = as.integer(barrier),
    spde_barrier = make_barrier_spde(spde),
    barrier_scaling = if (barrier) spde$barrier_scaling else c(1, 1),
    anisotropy = as.integer(anisotropy),
    family     = .valid_family[family$family],
    size = c(size),
    link       = .valid_link[family$link],
    df         = df,
    spatial_only = as.integer(spatial_only),
    spatial_covariate = as.integer(!is.null(spatial_varying)),
    calc_quadratic_range = as.integer(quadratic_roots),
    X_threshold = thresh$X_threshold,
    proj_X_threshold = 0, # dummy
    threshold_func = thresh$threshold_func,
    RE_indexes = RE_indexes,
    proj_RE_indexes = matrix(0, ncol = 0, nrow = 1), # dummy
    nobs_RE = nobs_RE,
    ln_tau_G_index = ln_tau_G_index,
    est_epsilon_model = as.integer(est_epsilon_model),
    epsilon_predictor = epsilon_covariate,
    est_epsilon_slope = as.integer(est_epsilon_slope),
    est_epsilon_re = as.integer(est_epsilon_re),
    has_smooths = as.integer(sm$has_smooths),
    upr = upr,
    lwr = lwr,
    poisson_link_delta = as.integer(isTRUE(family$type == "poisson_link_delta"))
  )

  b_thresh <- rep(0, 2)
  if (thresh$threshold_func == 2L) b_thresh <- c(0, b_thresh) # logistic

  tmb_params <- list(
    ln_H_input = c(0, 0),
    b_j        = matrix(0, ncol(X_ij), n_m),
    bs         = if (sm$has_smooths) matrix(0, nrow = ncol(sm$Xs), ncol = n_m) else array(0),
    ln_tau_O   = rep(0, n_m),
    ln_tau_Z = matrix(0, n_z, n_m),
    ln_tau_E   = rep(0, n_m),
    ln_kappa   = matrix(0, 2L, n_m),
    # ln_kappa   = rep(log(sqrt(8) / median(stats::dist(spde$mesh$loc))), 2),
    thetaf     = 0,
    ln_phi     = rep(0, n_m),
    ln_tau_V   = matrix(0, ncol(X_rw_ik), n_m),
    ar1_phi    = rep(0, n_m),
    ln_tau_G   = matrix(0, ncol(RE_indexes), n_m),
    RE         = matrix(0, sum(nobs_RE), n_m),
    b_rw_t     = array(0, dim = c(tmb_data$n_t, ncol(X_rw_ik), n_m)),
    omega_s    = matrix(0, n_s, n_m),
    zeta_s    = array(0, dim = c(n_s, n_z, n_m)),
    epsilon_st = array(0, dim = c(n_s, tmb_data$n_t, n_m)),
    b_threshold = b_thresh,
    # b_epsilon = 0,
    # ln_epsilon_re_sigma = 0,
    # epsilon_re = rep(0, tmb_data$n_t),
    b_smooth = if (sm$has_smooths) matrix(0, sum(sm$sm_dims), n_m) else array(0),
    ln_smooth_sigma = if (sm$has_smooths) matrix(0, length(sm$sm_dims), n_m) else array(0)
  )
  if (identical(family$link, "inverse") && family$family[1] %in% c("Gamma", "gaussian", "student") && !delta) {
    fam <- family
    if (family$family == "student") fam$family <- "gaussian"
    temp <- mgcv::gam(formula = formula, data = data, family = fam)
    tmb_params$b_j <- matrix(stats::coef(temp), ncol = 1L)
  }
  if (contains_offset) tmb_params$b_j[offset_pos] <- 1

  # Mapping off params as needed:
  tmb_map <- list()
  tmb_map$ln_phi <- rep(1, n_m)
  if (!anisotropy)
    tmb_map <- c(tmb_map, list(ln_H_input = factor(rep(NA, 2))))
  tmb_map$ar1_phi <- factor(rep(NA, n_m))
  if (family$family[[1]] %in% c("binomial", "poisson", "censored_poisson"))
    tmb_map$ln_phi[1] <- factor(NA)
  if (delta) {
    if (family$family[[2]] %in% c("binomial", "poisson", "censored_poisson"))
      tmb_map$ln_phi[2] <- factor(NA)
    else
      tmb_map$ln_phi[2] <- 2
  }
  tmb_map$ln_phi <- as.factor(tmb_map$ln_phi)

  if (!"tweedie" %in% family$family)
    tmb_map <- c(tmb_map, list(thetaf = as.factor(NA)))
  if (all(spatial_only))
    tmb_map <- c(tmb_map, list(
      ln_tau_E   = factor(rep(NA, n_m)),
      epsilon_st = factor(rep(NA, length(tmb_params$epsilon_st)))))

  if (delta && contains_offset) stop("`offset` in formula is deprecated in favour of the `offset` argument. For delta models, the `offset` argument must be used.", call. = FALSE)
  if (contains_offset) { # fix offset param to 1 to be an offset:
    b_j_map <- seq_along(tmb_params$b_j)
    b_j_map[offset_pos] <- NA
    tmb_map <- c(tmb_map, list(b_j = as.factor(b_j_map)))
  }

  if (delta && !is.null(thresh$threshold_parameter)) stop("Offsets not implemented with threshold models yet!") # TODO DELTA
  if (is.null(thresh$threshold_parameter)) {
    tmb_map <- c(tmb_map, list(b_threshold = factor(rep(NA, 2))))
  }

##  # optional models on spatiotemporal sd parameter
##  if (est_epsilon_re == 0L) {
##    tmb_map <- c(tmb_map, list(ln_epsilon_re_sigma = as.factor(NA),
##      epsilon_re = factor(rep(NA, tmb_data$n_t))))
##  }
##  if (est_epsilon_model == 0L) {
##    tmb_map <- c(tmb_map, list(b_epsilon = as.factor(NA)))
##  }

  if (multiphase && is.null(previous_fit) && do_fit) {
    not_phase1 <- c(tmb_map, list(
      ln_tau_O   = factor(rep(NA, n_m)),
      ln_tau_E   = factor(rep(NA, n_m)),
      ln_tau_V   = factor(rep(NA, length(tmb_params$ln_tau_V))),
      ln_tau_G   = factor(rep(NA, length(tmb_params$ln_tau_G))),
      ln_tau_Z   = factor(rep(NA, length(tmb_params$ln_tau_Z))),
      zeta_s     = factor(rep(NA, length(tmb_params$zeta_s))),
      ln_kappa   = factor(rep(NA, length(tmb_params$ln_kappa))),
      ln_H_input = factor(rep(NA, 2)),
      b_rw_t     = factor(rep(NA, length(tmb_params$b_rw_t))),
      RE         = factor(rep(NA, length(tmb_params$RE))),
      omega_s    = factor(rep(NA, length(tmb_params$omega_s))),
      b_smooth = factor(rep(NA, length(tmb_params$b_smooth))),
      ln_smooth_sigma = factor(rep(NA, length(tmb_params$ln_smooth_sigma))),
      epsilon_st = factor(rep(NA, length(tmb_params$epsilon_st)))))

    tmb_obj1 <- TMB::MakeADFun(
      data = tmb_data, parameters = tmb_params,
      map = not_phase1, DLL = "sdmTMB", silent = silent)

    lim <- set_limits(tmb_obj1, lower = lower, upper = upper, silent = TRUE)
    tmb_opt1 <- stats::nlminb(
      start = tmb_obj1$par, objective = tmb_obj1$fn,
      lower = lim$lower, upper = lim$upper,
      gradient = tmb_obj1$gr, control = .control)

##    # Set starting values based on phase 1:
##    if (isFALSE(contains_offset))
##      tmb_params$b_j <- set_par_value(tmb_opt1, "b_j") # TODO DELTA
##    else
##      tmb_params$b_j[-offset_pos] <- set_par_value(tmb_opt1, "b_j") # TODO DELTA

    if (family$family[[1]] == "tweedie")
      tmb_params$thetaf <- set_par_value(tmb_opt1, "thetaf")
    if (!family$family[[1]] %in% c("binomial", "poisson", "censored_poisson") && !delta)  # no dispersion param
      tmb_params$ln_phi <- set_par_value(tmb_opt1, "ln_phi")
  }

  if (all(spatial_only)) {
    tmb_random <- "omega_s"
  } else {
    if (include_spatial) {
      tmb_random <- c("omega_s", "epsilon_st")
    } else {
      tmb_random <- "epsilon_st"
    }
  }
  if (!is.null(spatial_varying)) tmb_random <- c(tmb_random, "zeta_s")
  if (!is.null(time_varying)) tmb_random <- c(tmb_random, "b_rw_t")
  if(est_epsilon_re) tmb_random <- c(tmb_random, "epsilon_re")

  if (!include_spatial) {
    tmb_map <- c(tmb_map, list(
      ln_tau_O = factor(rep(NA, n_m)),
      omega_s  = factor(rep(NA, length(tmb_params$omega_s)))))
  }
  if (is.null(spatial_varying)) {
    tmb_map <- c(tmb_map, list(
      ln_tau_Z = factor(rep(NA, n_m * n_z)),
      zeta_s = factor(rep(NA, length(tmb_params$zeta_s)))))
  }
  if (is.null(time_varying))
    tmb_map <- c(tmb_map,
      list(b_rw_t = factor(rep(NA, length(tmb_params$b_rw_t)))),
      list(ln_tau_V = factor(rep(NA, n_m)))
    )

  tmb_map$ar1_phi <- rep(NA, n_m)
  for (i in seq_along(spatiotemporal)) {
    if (spatiotemporal[i] == "ar1") tmb_map$ar1_phi[i] <- i
  }
  tmb_map$ar1_phi <- as.factor(as.integer(as.factor(tmb_map$ar1_phi)))

  if (nobs_RE[[1]] > 0) tmb_random <- c(tmb_random, "RE")
  if (reml) tmb_random <- c(tmb_random, "b_j")

##  if (est_epsilon_model >= 2) {
##    # model 2 = re model, model 3 = loglinear-re
##    tmb_random <- c(tmb_random, "epsilon_rw")
##  }

  if (sm$has_smooths) {
    tmb_random <- c(tmb_random, "b_smooth") # smooth random effects
    # warning(
    #   "Detected a `s()` smoother. Smoothers are penalized in sdmTMB as\n",
    #   "of version 0.0.19, but used to be unpenalized.\n",
    #   "You no longer need to specify `k` since the degree of wiggliness\n",
    #   "is determined by the data.", call. = FALSE)
  }

  if (!is.null(previous_fit)) {
    tmb_params <- previous_fit$tmb_obj$env$parList()
  }

  tmb_data$normalize_in_r <- as.integer(normalize)

  if (!is.null(previous_fit)) tmb_map <- previous_fit$tmb_map
  if (isTRUE(map_rf)) tmb_map <- map_off_rf(tmb_map, tmb_params)
  tmb_map <- c(map, tmb_map)

  # if (share_range && !delta) tmb_map <- c(tmb_map, list(ln_kappa = as.factor(rep(1, length(tmb_params$ln_kappa)))))
  # if (share_range && delta) tmb_map <- c(tmb_map, list(ln_kappa = as.factor(c(1, 1, 2, 2))))

  # deal with kappa mapping:
  tmb_map$ln_kappa <- matrix(seq_len(length(tmb_params$ln_kappa)),
    nrow(tmb_params$ln_kappa), ncol(tmb_params$ln_kappa))
  for (m in seq_len(n_m)) {
    if (share_range[m]) tmb_map$ln_kappa[,m] <- tmb_map$ln_kappa[1,m]
    if (spatiotemporal[m] == "off") tmb_map$ln_kappa[,m] <- NA
  }
  tmb_map$ln_kappa <- as.factor(as.integer(as.factor(tmb_map$ln_kappa)))

  for (i in seq_along(start)) {
    message("Initiating ", names(start)[i],
      " at specified starting value of ", start[[i]], ".")
    tmb_params[[names(start)[i]]] <- start[[i]]
  }

  if (nrow(tmb_params[["ln_kappa"]]) != 2L && "ln_kappa" %in% names(start)) {
    stop("Note that `ln_kappa` must be a matrix of nrow 2 and ncol models (regular=1, delta=2). It should be the same value in each row if `share_range = TRUE`.", call. = FALSE)
  }

  data$sdm_x <- data$sdm_y <- data$sdm_orig_id <- data$sdm_spatial_id <- NULL
  data$sdmTMB_X_ <- data$sdmTMB_Y_ <- NULL

  if (!tmb_data$has_smooths) {
    tmb_map <- c(tmb_map, list(b_smooth = factor(NA)))
    tmb_map <- c(tmb_map, list(bs = factor(NA)))
    tmb_map <- c(tmb_map, list(ln_smooth_sigma = factor(NA)))
  }

  # FIXME: generalize ; DELTA
  if (delta && "off" %in% spatiotemporal) {
    tmb_map$epsilon_st <- array(
      seq_len(length(tmb_params$epsilon_st)),
      dim = dim(tmb_params$epsilon_st)
    )
    tmb_map$ln_tau_E <- seq_len(length(tmb_params$ln_tau_E))
    for (i in which(spatiotemporal == "off")) {
      tmb_map$epsilon_st[,,i] <- NA
      tmb_map$ln_tau_E[i] <- NA
    }
    tmb_map$epsilon_st <- as.factor(tmb_map$epsilon_st)
    tmb_map$ln_tau_E <- as.factor(tmb_map$ln_tau_E)
  }

  tmb_obj <- TMB::MakeADFun(
    data = tmb_data, parameters = tmb_params, map = tmb_map,
    profile = if (control$profile) "b_j" else NULL,
    random = tmb_random, DLL = "sdmTMB", silent = silent)
  lim <- set_limits(tmb_obj, lower = lower, upper = upper,
    loc = spde$mesh$loc, silent = FALSE)

  out_structure <- structure(list(
    data       = data,
    spde       = spde,
    mesh       = spde,
    formula    = original_formula,
    split_formula = split_formula,
    time_varying = time_varying,
    threshold_parameter = thresh$threshold_parameter,
    threshold_function = thresh$threshold_func,
    epsilon_predictor = epsilon_predictor,
    time       = time,
    family     = family,
    response   = y_i,
    tmb_data   = tmb_data,
    tmb_params = tmb_params,
    tmb_obj    = tmb_obj,
    tmb_map    = tmb_map,
    tmb_random = tmb_random,
    spatial_varying = spatial_varying,
    spatial_varying_formula = spatial_varying_formula,
    reml       = reml,
    lower      = lim$lower,
    upper      = lim$upper,
    priors     = priors,
    nlminb_control = .control,
    control  = control,
    contrasts  = attr(X_ij, "contrasts"),
    terms  = attr(mf, "terms"),
    extra_time = extra_time,
    xlevels    = stats::.getXlevels(mt, mf),
    call       = match.call(expand.dots = TRUE),
    version    = utils::packageVersion("sdmTMB")),
    class      = "sdmTMB")
  if (!do_fit) return(out_structure)

  if (normalize) tmb_obj <- TMB::normalize(tmb_obj, flag = "flag", value = 0)

  tmb_opt <- stats::nlminb(
    start = tmb_obj$par, objective = tmb_obj$fn, gradient = tmb_obj$gr,
    lower = lim$lower, upper = lim$upper, control = .control)

  if (nlminb_loops > 1) {
    if(!silent) cat("running extra nlminb loops\n")
    for (i in seq(2, nlminb_loops, length = max(0, nlminb_loops - 1))) {
      temp <- tmb_opt[c("iterations", "evaluations")]
      tmb_opt <- stats::nlminb(
        start = tmb_opt$par, objective = tmb_obj$fn, gradient = tmb_obj$gr,
        control = .control, lower = lim$lower, upper = lim$upper)
      tmb_opt[["iterations"]] <- tmb_opt[["iterations"]] + temp[["iterations"]]
      tmb_opt[["evaluations"]] <- tmb_opt[["evaluations"]] + temp[["evaluations"]]
    }
  }
  if (newton_loops > 0) {
    if(!silent) cat("running newtonsteps\n")
    for (i in seq_len(newton_loops)) {
      g <- as.numeric(tmb_obj$gr(tmb_opt$par))
      h <- stats::optimHess(tmb_opt$par, fn = tmb_obj$fn, gr = tmb_obj$gr)
      tmb_opt$par <- tmb_opt$par - solve(h, g)
      tmb_opt$objective <- tmb_obj$fn(tmb_opt$par)
    }
  }
  check_bounds(tmb_opt$par, lim$lower, lim$upper)

  sd_report <- TMB::sdreport(tmb_obj, getJointPrecision = get_joint_precision)
  conv <- get_convergence_diagnostics(sd_report)

  out_structure$tmb_obj <- tmb_obj
  out <- c(out_structure, list(
    model      = tmb_opt,
    sd_report  = sd_report,
    gradients  = conv$final_grads,
    bad_eig    = conv$bad_eig,
    pos_def_hessian = sd_report$pdHess))
  `class<-`(out, "sdmTMB")
}

map_off_rf <- function(.map, tmb_params) {
  .map$ln_tau_O <- as.factor(rep(NA, length(tmb_params$ln_tau_O)))
  .map$ln_tau_E <- as.factor(rep(NA, length(tmb_params$ln_tau_E)))
  .map$ln_tau_Z <- as.factor(rep(NA, length(tmb_params$ln_tau_Z)))
  .map$zeta_s <- factor(rep(NA, length(tmb_params$zeta_s)))
  .map$ln_kappa <- factor(rep(NA, length(tmb_params$ln_kappa)))
  .map$ln_H_input <- factor(rep(NA, length(tmb_params$ln_H_input)))
  .map$omega_s <- factor(rep(NA, length(tmb_params$omega_s)))
  .map$epsilon_st <- factor(rep(NA, length(tmb_params$epsilon_st)))
  .map
}

check_bounds <- function(.par, lower, upper) {
  for (i in seq_along(.par)) {
    if (is.finite(lower[i]) || is.finite(upper[i])) {
      if (abs(.par[i] - lower[i]) < 0.0001) {
        warning("Parameter ", names(.par)[i],
          " is very close or equal to its lower bound.\n",
          "Consider changing your model configuration or bounds.",
          call. = FALSE)
      }
      if (abs(.par[i] - upper[i]) < 0.0001) {
        warning("Parameter ", names(.par)[i],
          " is very close or equal to its upper bound.\n",
          "Consider changing your model configuration or bounds.",
          call. = FALSE)
      }
    }
  }
}

set_limits <- function(tmb_obj, lower, upper, loc = NULL, silent = TRUE) {
  .lower <- stats::setNames(rep(-Inf, length(tmb_obj$par)), names(tmb_obj$par))
  .upper <- stats::setNames(rep(Inf, length(tmb_obj$par)), names(tmb_obj$par))
  for (i_name in names(lower)) {
    if (i_name %in% names(.lower)) {
      .lower[names(.lower) %in% i_name] <- lower[[i_name]]
      if (!silent) message("Setting lower limit for ", i_name, " to ",
        lower[[i_name]], ".")
    }
    if (i_name %in% names(.upper)) {
      .upper[names(.upper) %in% i_name] <- upper[[i_name]]
      if (!silent) message("Setting upper limit for ", i_name, " to ",
        upper[[i_name]], ".")
    }
  }
  if ("ar1_phi" %in% names(tmb_obj$par) &&
      !"ar1_phi" %in% union(names(lower), names(upper))) {
    .lower["ar1_phi"] <- stats::qlogis((-0.999 + 1) / 2)
    .upper["ar1_phi"] <- stats::qlogis((0.999 + 1) / 2)
  }

  # if (!is.null(loc) && !"ln_kappa" %in% union(names(lower), names(upper))) {
  #   .dist <- stats::dist(loc)
  #   range_limits <- c(min(.dist) * 0.25, max(.dist) * 4)
  #   .upper[names(.upper) == "ln_kappa"] <- log(sqrt(8) / range_limits[1])
  #   .lower[names(.upper) == "ln_kappa"] <- log(sqrt(8) / range_limits[2])
  #   message("Setting limits on range to ",
  #     round(range_limits[1], 1), " and ",
  #     round(range_limits[2], 1), ":\n",
  #     "half the minimum and twice the maximum knot distance."
  #   )
  #   if (.upper[names(.upper) == "ln_kappa"][1] <= tmb_obj$par[["ln_kappa"]][1]) {
  #     .upper[names(.upper) == "ln_kappa"] <- tmb_obj$par[["ln_kappa"]][1] + 0.1
  #   }
  #   if (.lower[names(.lower) == "ln_kappa"][1] >= tmb_obj$par[["ln_kappa"]][1]) {
  #     .lower[names(.lower) == "ln_kappa"] <- tmb_obj$par[["ln_kappa"]][1] - 0.1
  #   }
  # }
  list(lower = .lower, upper = .upper)
}

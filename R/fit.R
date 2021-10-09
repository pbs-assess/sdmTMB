#' @useDynLib sdmTMB, .registration = TRUE
NULL

#' Fit a spatial or spatiotemporal GLMM with TMB
#'
#' Fit a spatial or spatiotemporal predictive-process GLMM with TMB. This can be
#' useful for (dynamic) species distribution models and relative abundance index
#' standardization among many other uses.
#'
#' @param formula Model formula. See the Details section below for how to specify
#'   offsets and threshold parameters. For index standardization, you may wish
#'   to include `0 + as.factor(year)` (or whatever the time column is called)
#'   in the formula. IID random intercepts are possible using \pkg{lme4}
#'   syntax, e.g., `+ (1 | g)` where `g` is a column with factor levels.
#'   Penalized splines are possible via \pkg{mgcv}. See examples below.
#' @param data A data frame.
#' @param spde An object from [make_mesh()].
#' @param time An optional time column name (as character). Can be left as
#'   `NULL` for a model with only spatial random fields unless you wish to use
#'   one of the index or center of gravity functions over time.
#' @param family The family and link. Supports [gaussian()], [Gamma()],
#'   [binomial()], [poisson()], \code{\link[sdmTMB:families]{Beta()}},
#'   \code{\link[sdmTMB:families]{nbinom2()}}, and
#'   \code{\link[sdmTMB:families]{tweedie()}}]. For binomial family options,
#'   see the 'Binomial families' in the Details section below.
#' @param time_varying An optional formula describing covariates that should be
#'   modelled as a random walk through time. Be careful not to include
#'   covariates (including the intercept) in both the main and time-varying
#'   formula. I.e., at least one should have `~ 0` or `~ -1`.
#' @param weights Optional likelihood weights for the conditional model.
#'   Implemented as in \pkg{glmmTMB}. In other words, weights do not have to sum
#'   to one and are not internally modified. Can also be used for trials with
#'   the binomial family. See Details below.
#' @param extra_time Optional extra time slices (e.g., years) to include for
#'   interpolation or forecasting with the predict function. See the
#'   Details section below.
#' @param reml Logical: use REML (restricted maximum likelihood) estimation
#'   rather than maximum likelihood? Internally, this adds the fixed effects
#'   to the list of random effects to intetgrate over.
#' @param silent Silent or include optimization details?
#' @param anisotropy Logical: allow for anisotropy? See [plot_anisotropy()].
#' @param control Optimization control options. See [sdmTMBcontrol()].
#' @param priors Optional penalties/priors. See [sdmTMBpriors()].
#' @param fields Estimate the spatiotemporal random fields as `"IID"`
#'   (independent and identically distributed; default), stationary `"AR1"`
#'   (first-order autoregressive), or as a random walk (`"RW"`). Note that the
#'   spatiotemporal standard deviation represents the marginal steady-state
#'   standard deviation of the process in the case of the AR1. I.e., it is
#'   scaled according to the correlation. See the [TMB
#'   documentation](https://kaskr.github.io/adcomp/classAR1__t.html). If the AR1
#'   correlation coefficient (rho) is estimated close to 1, say > 0.99, then you
#'   may wish to switch to the random walk `"RW"`.
#' @param include_spatial Should a separate spatial random field be estimated?
#'   If enabled then there will be separate spatial and spatiotemporal
#'   fields.
#' @param spatial_trend Should a separate spatial field be included in the
#'   trend that represents local time trends? Requires spatiotemporal data.
#'   See \doi{10.1111/ecog.05176} and the
#'   [spatial trends vignette](https://pbs-assess.github.io/sdmTMB/articles/spatial-trend-models.html).
#' @param spatial_only Logical: should only a spatial model be fit (i.e., do not
#'   include spatiotemporal random effects)? By default a spatial-only model
#'   will be fit if there is only one unique value in the time column or the
#'   `time` argument is left at its default value of `NULL`.
#' @param share_range Logical: estimate a shared spatial and spatiotemporal
#'   range parameter (`TRUE`) or independent range parameters (`FALSE`).
#' @param previous_fit A previously fitted sdmTMB model to initialize the
#'   optimization with. Can greatly speed up fitting. Note that the data and
#'   model must be set up *exactly* the same way. However, the `weights` argument
#'   can change, which can be useful for cross-validation.
#' @param epsilon_predictor (Experimental) A column name (as character) of a
#'   predictor of a linear trend (in log space) of the spatiotemporal standard
#'   deviation. By default, this is `NULL` and fits a model with a constant
#'   spatiotemporal variance. However, this argument can also be a character
#'   name in the original data frame (a covariate that ideally has been
#'   standardized to have mean 0 and standard deviation = 1). Because the
#'   spatiotemporal field varies by time step, the standardization should be
#'   done by time. If the name of a predictor is included, a log-linear model is
#'   fit where the predictor is used to model effects on the standard deviation,
#'   e.g. `log(sd(i)) = B0 + B1 * epsilon_predictor(i)`.
#' @param do_fit Fit the model (`TRUE`) or return the processed data without
#'   fitting (`FALSE`)?
#' @param ... Not currently used.
#' @importFrom methods as is
#' @importFrom mgcv s t2
#' @importFrom stats gaussian model.frame model.matrix
#'   model.response terms model.offset
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
#' **Offsets**
#'
#' In the model formula, an offset can be included by including `+ offset` in
#' the model formula (a reserved word). The offset will be included in any
#' prediction. `offset` must be a column in `data`.
#'
#' **Binomial families**

#' Following the structure of [stats::glm()] and \pkg{glmmTMB}, a binomial
#' family can be specified in one of 4 ways: (1) the response may be a factor
#' (and the model classifies the first level versus all others), (2) the
#' response may be binomial (0/1), (3) the response can be a matrix of form
#' `cbind(success, failure)`, and (4) the response may be the observed
#' proportions, and the 'weights' argument is used to specify the Binomial size
#' (N) parameter (`prob ~ ..., weights = N`).
#'
#' **Threshold models**
#'
#' A linear break-point relationship for a covariate can be included via
#' `+ breakpt(variable)` in the formula, where `variable` is a single covariate
#' corresponding to a column in `data`. In this case, the relationship is linear
#' up to a point and then constant.
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
#' **Forecasting or interpolating**
#'
#' Extra time slices (e.g., years) can be included for interpolation or
#' forecasting with the predict function via the `extra_time` argument. The
#' predict function requires all time slices to be defined when fitting the
#' model to ensure the various time indices are set up correctly. Be careful if
#' including extra time slices that the model remains identifiable. For example,
#' including `+ as.factor(year)` in `formula` will render a model with no data
#' to inform the expected value in a missing year. [sdmTMB()] makes no attempt
#' to determine if the model makes sense for forecasting or interpolation. The
#' options `time_varying`, `fields = "RW"`, and `fields = "AR1"`
#' provide mechanisms to predict over missing time slices with process error.
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
#' parameters. See [sdmTMBpriors()]. These should not include `offset` terms and
#' care should be taken if used with splines. You can fit the model once without
#' penalties and inspect the element `head(your_model$tmb_data$X_ij)` if you
#' want to see how the formula is translated to the fixed effect model matrix.
#'
#' @references
#'
#' Main reference/report introducing the package. We plan to write a paper
#' to cite in the near future:
#'
#' Anderson, S.C., E.A. Keppel, A.M. Edwards, 2019. A reproducible data synopsis
#' for over 100 species of British Columbia groundfish. DFO Can. Sci. Advis. Sec.
#' Res. Doc. 2019/041. vii + 321 p.
#' <https://www.dfo-mpo.gc.ca/csas-sccs/Publications/ResDocs-DocRech/2019/2019_041-eng.html>
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
#' A.M. Edwards, B.M. Connors, S.C. Anderson. Contrasting climate velocity impacts
#' in warm and cool locations: A meta-analysis across 38 demersal fish species in
#' the northeast Pacific. EcoEvoRxiv. \doi{10.32942/osf.io/b87ng}.
#'
#' Code for implementing the penalized complexity prior on the spatial
#' covariance parameters:
#'
#' Osgood-Zimmerman, A. and Wakefield, J. 2021. A Statistical introduction
#' to Template Model Builder: a flexible tool for spatial modeling.
#' arXiv 2103.09929. <https://arxiv.org/abs/2103.09929>.
#'
#' Code for implementing the barrier-SPDE written by Olav Nikolai Breivik and
#' Hans Skaug and adapted via the VAST R package.
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
#'   pcod_spde <- make_mesh(pcod_2011, c("X", "Y"), cutoff = 25) # a coarse mesh for example speed
#' } else {
#'   pcod_spde <- pcod_mesh_2011 # load cached mesh so example runs on CRAN
#' }
#' plot(pcod_spde)
#'
#' # Tweedie:
#' m <- sdmTMB(density ~ 0 + depth_scaled + depth_scaled2 + as.factor(year),
#'   data = pcod_2011, time = "year", spde = pcod_spde, family = tweedie(link = "log"))
#' print(m)
#' tidy(m, conf.int = TRUE)
#' tidy(m, effects = "ran_par", conf.int = TRUE)
#'
#' if (inla_installed()) {
#' # Run extra optimization steps to help convergence:
#' # (Not typically needed)
#' m1 <- run_extra_optimization(m, nlminb_loops = 0, newton_loops = 1)
#' max(m$gradients)
#' max(m1$gradients)
#'
#' # Bernoulli:
#' pcod_binom <- pcod_2011
#' pcod_binom$present <- ifelse(pcod_binom$density > 0, 1L, 0L)
#' m_bin <- sdmTMB(present ~ 0 + as.factor(year) + depth_scaled + depth_scaled2,
#'   data = pcod_binom, time = "year", spde = pcod_spde,
#'   family = binomial(link = "logit"))
#' print(m_bin)
#'
#' # Fit a spatial-only model (by not specifying `time`):
#' m <- sdmTMB(
#'   density ~ depth_scaled + depth_scaled2, data = pcod_2011,
#'   spde = pcod_spde, family = tweedie(link = "log"))
#' print(m)
#'
#' # Gaussian:
#' pcod_gaus <- subset(pcod_2011, density > 0 & year >= 2013)
#' pcod_spde_gaus <- make_mesh(pcod_gaus, c("X", "Y"), cutoff = 30)
#' m_pos <- sdmTMB(log(density) ~ 0 + as.factor(year) + depth_scaled + depth_scaled2,
#'   data = pcod_gaus, time = "year", spde = pcod_spde_gaus)
#' print(m_pos)
#'
#' # With p-splines via mgcv:
#' m_gam <- sdmTMB(log(density) ~ s(depth_scaled),
#'   data = pcod_gaus, time = "year", spde = pcod_spde_gaus)
#'
#' # With IID random intercepts:
#' # Simulate some data:
#' set.seed(1)
#' x <- runif(500, -1, 1)
#' y <- runif(500, -1, 1)
#' loc <- data.frame(x = x, y = y)
#' spde <- make_mesh(loc, c("x", "y"), n_knots = 50, type = "kmeans")
#' s <- sdmTMB_sim(x = x, y = y, betas = 0, time = 1L,
#'   phi = 0.1, range = 1.4, sigma_O = 0.2, sigma_E = 0, mesh = spde)
#' s$g <- gl(50, 10)
#' iid_re_vals <- rnorm(50, 0, 0.3)
#' s$observed <- s$observed + iid_re_vals[s$g]
#'
#' # Fit it:
#' m <- sdmTMB(observed ~ 1 + (1 | g), spde = spde, data = s)
#' print(m)
#' tidy(m, "ran_pars", conf.int = TRUE) # see tau_G
#' theta <- as.list(m$sd_report, "Estimate")
#' plot(iid_re_vals, theta$RE)
#' }
#'
#' \donttest{
#' if (inla_installed()) {
#'
#' # Spatial-trend example:
#' m <- sdmTMB(density ~ depth_scaled, data = pcod_2011,
#'   spde = pcod_spde, family = tweedie(link = "log"),
#'   spatial_trend = TRUE, time = "year")
#' tidy(m, effects = "ran_par")
#'
#' # Time-varying effects of depth and depth squared:
#' m <- sdmTMB(density ~ 0 + as.factor(year),
#'   time_varying = ~ 0 + depth_scaled + depth_scaled2,
#'   data = pcod_2011, time = "year", spde = pcod_spde,
#'   family = tweedie(link = "log"))
#' print(m)
#'
#' # See the b_rw_t estimates; these are the time-varying (random walk) effects.
#' # These could be added to tidy.sdmTMB() eventually.
#' summary(m$sd_report)[1:19,]
#'
#' # Linear breakpoint model on depth:
#' m_pos <- sdmTMB(log(density) ~ 0 + as.factor(year) +
#'     breakpt(depth_scaled) + depth_scaled2, data = pcod_gaus,
#'   time = "year", spde = pcod_spde_gaus)
#' print(m_pos)
#'
#' # Linear covariate on log(sigma_epsilon):
#' # First we will center the years around their mean
#' # to help with convergence. Also constrain the slope to be (-1, 1)
#' # to help convergence (estimation is done in log-space)
#' pcod_2011$year_centered <- pcod_2011$year - mean(pcod_2011$year)
#' m <- sdmTMB(density ~ 0 + depth_scaled + depth_scaled2 + as.factor(year),
#'   data = pcod_2011, time = "year", spde = pcod_spde, family = tweedie(link = "log"),
#'   epsilon_predictor = "year_centered",
#'   control = sdmTMBcontrol(lower = list(b_epsilon = -1), upper = list(b_epsilon = 1)))
#' print(m) # sigma_E varies with time now
#'
#' # coefficient is not yet in tidy.sdmTMB:
#' as.list(m$sd_report, "Estimate", report = TRUE)$b_epsilon
#' as.list(m$sd_report, "Std. Error", report = TRUE)$b_epsilon
#' }
#' }

sdmTMB <- function(formula, data, spde, time = NULL,
  family = gaussian(link = "identity"),
  time_varying = NULL, weights = NULL, extra_time = NULL, reml = FALSE,
  silent = TRUE, anisotropy = FALSE,
  control = sdmTMBcontrol(),
  priors = sdmTMBpriors(),
  fields = c("IID", "AR1", "RW"),
  include_spatial = TRUE,
  spatial_trend = FALSE,
  spatial_only = identical(length(unique(data[[time]])), 1L),
  share_range = TRUE,
  previous_fit = NULL,
  epsilon_predictor = NULL,
  do_fit = TRUE,
  ...
  ) {

  normalize <- control$normalize
  nlminb_loops <- control$nlminb_loops
  newton_loops <- control$newton_loops
  mgcv <- control$mgcv
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
    warning("`ar1_fields` is depreciated and is now specified via the ",
      "`fields` argument. Setting `fields = 'AR1'` for you for now if `ar1_fields = TRUE`.",
      call. = FALSE)
    fields <- if (dots$ar1_fields) "AR1" else "IID"
    dots$ar1_fields <- NULL
  }
  if ("penalties" %in% names(dots)) {
    stop("`penalties` are now specified via the `priors` argument.",
      "E.g., `priors = sdmTMBpriors(b = normal(c(0, 0), c(1, 1)))`",
      "for 2 fixed effects.", call. = FALSE)
  }
  dot_checks <- c("lower", "upper", "profile",
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

  fields <- match.arg(fields)
  ar1_fields <- identical("AR1", fields)
  rw_fields <- identical("RW", fields)
  assert_that(
    is.logical(reml), is.logical(anisotropy), is.logical(silent),
    is.logical(silent), is.logical(spatial_trend), is.logical(mgcv),
    is.logical(multiphase), is.logical(ar1_fields),
    is.logical(include_spatial), is.logical(map_rf), is.logical(normalize)
  )
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

  thresh <- check_and_parse_thresh_params(formula, data)
  original_formula <- formula
  formula <- thresh$formula

  if (!is.null(extra_time)) { # for forecasting or interpolating
    if (!"xy_cols" %in% names(spde)) {
      stop("Please use make_mesh() instead of the depreciated make_spde() to use `extra_time`.",
        call. = FALSE)
    }
    data <- expand_time(df = data, time_slices = extra_time, time_column = time)
    weights <- data$weight_sdmTMB
    spde$A <- INLA::inla.spde.make.A(spde$mesh, loc = as.matrix(data[, spde$xy_cols, drop = FALSE]))
    spde$loc_xy <- as.matrix(data[,spde$xy_cols,drop=FALSE])
  }

  if (spatial_trend) {
    numeric_time <- time
    t_i <- as.numeric(data[[numeric_time]])
    t_i <- t_i - mean(unique(t_i), na.rm = TRUE) # middle year = intercept
  } else {
    t_i <- rep(0L, nrow(data))
  }
  contains_offset <- check_offset(formula)

  split_formula <- glmmTMB::splitForm(formula)
  RE_names <- barnames(split_formula$reTrmFormulas)
  fct_check <- vapply(RE_names, function(x) check_valid_factor_levels(data[[x]], .name = x), TRUE)
  RE_indexes <- vapply(RE_names, function(x) as.integer(data[[x]]) - 1L, rep(1L, nrow(data)))
  nobs_RE <- unname(apply(RE_indexes, 2L, max)) + 1L
  if (length(nobs_RE) == 0L) nobs_RE <- 0L
  formula <- split_formula$fixedFormula
  ln_tau_G_index <- unlist(lapply(seq_along(nobs_RE), function(i) rep(i, each = nobs_RE[i]))) - 1L

  if (isFALSE(mgcv)) {
    mgcv_mod <- NULL
    X_ij <- model.matrix(formula, data)
    mf <- model.frame(formula, data)
  } else {
    # mgcv::gam will parse a matrix response, but not a factor
    mf <- model.frame(mgcv::interpret.gam(formula)$fake.formula, data)
    if (identical(family$family, "binomial") && "factor" %in% model.response(mf, "any")) {
      stop("Error: with 'mgcv' = TRUE, the response cannot be a factor", call. = FALSE)
    }
    if (family$family %in% c("binomial", "Gamma")) {
      mgcv_mod <- mgcv::gam(formula, data = data, family = family, fit = FALSE) # family needs to be passed into mgcv
    } else {
      mgcv_mod <- mgcv::gam(formula, data = data, fit = FALSE)
    }
    # X_ij <- model.matrix(mgcv_mod)
    X_ij <- mgcv_mod$X
    X_ij <- X_ij[, colnames(X_ij) != "", drop = FALSE] # only keep non-smooth terms
  }

  terms <- all_terms(formula)
  smooth_i <- get_smooth_terms(terms)
  basis <- list()
  rasm <- list()
  Zs <- list()
  Xs <- list()
  if (length(smooth_i) > 0) {
    has_smooths <- TRUE
    smterms <- terms[smooth_i]
    ns <- 0
    for (i in seq_along(smterms)) {
      obj <- eval(str2expression(smterms[i]))
      basis[[i]] <- mgcv::smoothCon(
        object = obj, data = data,
        knots = NULL, absorb.cons = TRUE, # modCon = 3, # see brms standata_basis_sm
        diagonal.penalty = TRUE
      )
      for (j in seq_along(basis[[i]])) { # basis is length > 1 with s(..., by = ...) terms
        ns <- ns + 1
        rasm[[ns]] <- mgcv::smooth2random(basis[[i]][[j]], names(data), type = 2)
        Zs[[ns]] <- rasm[[ns]]$rand$Xr
        Xs[[ns]] <- rasm[[ns]]$Xf
      }
    }
    sm_dims <- unlist(lapply(Zs, ncol))
    Xs <- do.call(cbind, Xs) # combine 'em all into one design matrix
    b_smooth_start <- c(0, cumsum(sm_dims)[-length(sm_dims)])
  } else {
    has_smooths <- FALSE
    sm_dims <- 0L
    b_smooth_start <- 0L
  }

  # FIXME?: deal with "by =" stuff
  # EW: thinking of the same approach brms takes with factors?
  # DONE: split off non-smooth model matrix
  # DONE: pass in Zs to TMB
  # DONE: set up weights coefs in R then TMB
  # DONE: add Normal() penalty on Zs weights in TMB
  # DONE: take normal deviations and multiply each sub-vector by corresponding sparse matrix
  # FIXME: deal with prediction (I think this is going to need additional matrices)
  # FIXME: test with s, t2, by, cc

  offset_pos <- grep("^offset$", colnames(X_ij))
  y_i <- model.response(mf, "numeric")
  if (family$family %in% c("Gamma", "lognormal") && min(y_i) <= 0) {
    stop("Gamma and lognormal must have response values > 0.", call. = FALSE)
  }

  # This is taken from approach in glmmTMB to match how they handle binomial
  # yobs could be a factor -> treat as binary following glm
  # yobs could be cbind(success, failure)
  # yobs could be binary
  # (yobs, weights) could be (proportions, size)
  # On the C++ side 'yobs' must be the number of successes.
  size <- rep(1, nrow(X_ij)) # for non-binomial case
  if (identical(family$family, "binomial")) {
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

  # if (!is.null(penalties)) {
  #   assert_that(ncol(X_ij) == length(penalties),
  #     msg = paste0("The number of fixed effects does not match the number of ",
  #       "penalty terms. Ensure that offset terms are not in penalty vector and ",
  #       "that any spline terms are properly accounted for."))
  # }

  if (identical(family$link, "log") && min(y_i, na.rm = TRUE) < 0) {
    stop("`link = 'log'` but the reponse data include values < 0.", call. = FALSE)
  }

  offset <- as.vector(model.offset(mf))
  if (is.null(offset)) offset <- rep(0, length(y_i))

  if (!is.null(time_varying)) {
    X_rw_ik <- model.matrix(time_varying, data)
  } else {
    X_rw_ik <- matrix(0, nrow = nrow(data), ncol = 1)
  }

  # Stuff needed for spatiotemporal A matrix:
  data$sdm_orig_id <- seq(1, nrow(data))
  data$sdm_x <- spde$loc_xy[,1,drop=TRUE]
  data$sdm_y <- spde$loc_xy[,2,drop=TRUE]
  fake_data <- unique(data.frame(sdm_x = data$sdm_x, sdm_y = data$sdm_y))
  fake_data[["sdm_spatial_id"]] <- seq(1, nrow(fake_data))
  data <- base::merge(data, fake_data, by = c("sdm_x", "sdm_y"),
    all.x = TRUE, all.y = FALSE)
  data <- data[order(data$sdm_orig_id),, drop=FALSE]
  A_st <- INLA::inla.spde.make.A(spde$mesh,
    loc = as.matrix(fake_data[, c("sdm_x", "sdm_y"), drop = FALSE]))

  n_s <- nrow(spde$mesh$loc)

  barrier <- "spde_barrier" %in% names(spde)
  if (barrier && anisotropy) {
    warning("Using a barrier mesh; therefore, anistropy will be disabled.", call. = FALSE)
    anisotropy <- FALSE
  }
  df <- if (family$family == "student" && "df" %in% names(family)) family$df else 3

  est_epsilon_model <- 0L
  epsilon_covariate <- rep(0, length(unique(data[[time]])))
  if (!is.null(epsilon_predictor)) {
    # covariate vector dimensioned by number of time steps
    time_steps <- unique(data[[time]])
    for (i in seq_along(time_steps)) {
      epsilon_covariate[i] <- data[data[[time]] == time_steps[i],
        epsilon_predictor, drop = TRUE][[1]]
    }
    est_epsilon_model <- 1L
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
    priors_b <- mvnormal(location = priors_b[, 1], scale = diag(as.numeric(priors_b[, 2])))
  }
  # in some cases, priors_b will be a mix of NAs (no prior) and numeric values
  # easiest way to deal with this is to subset the Sigma matrix
  not_na <- which(!is.na(priors_b[,1]))
  if (length(not_na) == 0L) {
    priors_b_Sigma <- diag(2) # TMB can't handle 0-dim matrix
  } else {
    Sigma <- priors_b[, -1]
    priors_b_Sigma <- as.matrix(Sigma[not_na, not_na])
  }

  tmb_data <- list(
    y_i        = c(y_i),
    n_t        = length(unique(data[[time]])),
    t_i        = t_i,
    offset_i   = offset,
    A          = spde$A,
    A_st       = A_st,
    A_spatial_index = data$sdm_spatial_id - 1L,
    year_i     = make_year_i(data[[time]]),
    ar1_fields = if (spatial_only) 0L else as.integer(ar1_fields),
    rw_fields = if (spatial_only) 0L else as.integer(rw_fields),
    X_ij       = X_ij,
    X_rw_ik    = X_rw_ik,
    Zs         = Zs, # optional smoother basis function matrices
    Xs         = Xs, # optional smoother linear effect matrix
    b_smooth_start = b_smooth_start,
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
    calc_time_totals = 0L,
    random_walk = !is.null(time_varying),
    priors_b_n = length(not_na),
    priors_b_index = not_na,
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
    proj_t_i = 0, # dummy
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
    spatial_trend = as.integer(spatial_trend),
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
    has_smooths = as.integer(has_smooths)
  )

  b_thresh <- rep(0, 2)
  if (thresh$threshold_func == 2L) b_thresh <- c(0, b_thresh) # logistic

  tmb_params <- list(
    ln_H_input = c(0, 0),
    b_j        = rep(0, ncol(X_ij)),
    bs         = rep(0, if (has_smooths) ncol(Xs) else 0),
    ln_tau_O   = 0,
    ln_tau_O_trend = 0,
    ln_tau_E   = 0,
    ln_kappa   = rep(0, 2),
    # ln_kappa   = rep(log(sqrt(8) / median(stats::dist(spde$mesh$loc))), 2),
    thetaf     = 0,
    ln_phi     = 0,
    ln_tau_V   = rep(0, ncol(X_rw_ik)),
    ar1_phi    = 0,
    ln_tau_G   = rep(0, ncol(RE_indexes)),
    RE         = rep(0, sum(nobs_RE)),
    b_rw_t     = matrix(0, nrow = tmb_data$n_t, ncol = ncol(X_rw_ik)),
    omega_s    = rep(0, n_s),
    omega_s_trend = rep(0, n_s),
    epsilon_st = matrix(0, nrow = n_s, ncol = tmb_data$n_t),
    b_threshold = b_thresh,
    b_epsilon = 0,
    b_smooth = rep(0, if (has_smooths) sum(sm_dims) else 0),
    ln_smooth_sigma = rep(0, if (has_smooths) length(sm_dims) else 0)
  )
  if (identical(family$link, "inverse") && family$family %in% c("Gamma", "gaussian", "student")) {
    fam <- family
    if (family$family == "student") fam$family <- "gaussian"
    temp <- mgcv::gam(formula = formula, data = data, family = fam)
    tmb_params$b_j <- stats::coef(temp)
  }
  if (contains_offset) tmb_params$b_j[offset_pos] <- 1

  # Mapping off params as needed:
  tmb_map <- list()
  if (!anisotropy)
    tmb_map <- c(tmb_map, list(ln_H_input = factor(rep(NA, 2))))
  if (!ar1_fields)
    tmb_map <- c(tmb_map, list(ar1_phi = as.factor(NA)))
  if (family$family %in% c("binomial", "poisson"))
    tmb_map <- c(tmb_map, list(ln_phi = as.factor(NA)))
  if (family$family != "tweedie")
    tmb_map <- c(tmb_map, list(thetaf = as.factor(NA)))
  if (spatial_only)
    tmb_map <- c(tmb_map, list(
      ln_tau_E   = as.factor(NA),
      epsilon_st = factor(rep(NA, length(tmb_params$epsilon_st)))))

  if (contains_offset) { # fix offset param to 1 to be an offset:
    b_j_map <- seq_along(tmb_params$b_j)
    b_j_map[offset_pos] <- NA
    tmb_map <- c(tmb_map, list(b_j = as.factor(b_j_map)))
  }

  if (is.null(thresh$threshold_parameter)) {
    tmb_map <- c(tmb_map, list(b_threshold = factor(rep(NA, 2))))
  }
  if (multiphase && is.null(previous_fit) && do_fit) {
    not_phase1 <- c(tmb_map, list(
      ln_tau_O   = as.factor(NA),
      ln_tau_E   = as.factor(NA),
      ln_tau_V   = factor(rep(NA, ncol(X_rw_ik))),
      ln_tau_G   = factor(rep(NA, length(tmb_params$ln_tau_G))),
      ln_tau_O_trend = as.factor(NA),
      omega_s_trend  = factor(rep(NA, length(tmb_params$omega_s_trend))),
      ln_kappa   = factor(rep(NA, 2)),
      ln_H_input = factor(rep(NA, 2)),
      b_rw_t     = factor(rep(NA, length(tmb_params$b_rw_t))),
      RE         = factor(rep(NA, length(tmb_params$RE))),
      omega_s    = factor(rep(NA, length(tmb_params$omega_s))),
      b_smooth = factor(rep(NA, length(tmb_params$b_smooth))),
      ln_smooth_sigma = factor(rep(NA, length(tmb_params$ln_smooth_sigma))),
      epsilon_st = factor(rep(NA, length(tmb_params$epsilon_st)))))

    # optional models on spatiotemporal sd parameter
    if (est_epsilon_model == 0L) {
      tmb_map <- c(tmb_map, list(b_epsilon = as.factor(NA)))
    }

    tmb_obj1 <- TMB::MakeADFun(
      data = tmb_data, parameters = tmb_params,
      map = not_phase1, DLL = "sdmTMB", silent = silent)

    lim <- set_limits(tmb_obj1, lower = lower, upper = upper, silent = TRUE)
    tmb_opt1 <- stats::nlminb(
      start = tmb_obj1$par, objective = tmb_obj1$fn,
      lower = lim$lower, upper = lim$upper,
      gradient = tmb_obj1$gr, control = .control)

    # Set starting values based on phase 1:
    if (isFALSE(contains_offset))
      tmb_params$b_j <- set_par_value(tmb_opt1, "b_j")
    else
      tmb_params$b_j[-offset_pos] <- set_par_value(tmb_opt1, "b_j")

    if (family$family == "tweedie")
      tmb_params$thetaf <- set_par_value(tmb_opt1, "thetaf")
    if (!family$family %in% c("binomial", "poisson"))  # no dispersion param
      tmb_params$ln_phi <- set_par_value(tmb_opt1, "ln_phi")
  }

  if (spatial_only) {
    tmb_random <- "omega_s"
  } else {
    if (include_spatial) {
      tmb_random <- c("omega_s", "epsilon_st")
    } else {
      tmb_random <- "epsilon_st"
    }
  }
  if (spatial_trend) tmb_random <- c(tmb_random, "omega_s_trend")
  if (!is.null(time_varying)) tmb_random <- c(tmb_random, "b_rw_t")

  if (!include_spatial) {
    tmb_map <- c(tmb_map, list(
      ln_tau_O = as.factor(NA),
      omega_s  = factor(rep(NA, length(tmb_params$omega_s)))))
  }
  if (!spatial_trend) {
    tmb_map <- c(tmb_map, list(
      ln_tau_O_trend = as.factor(NA),
      omega_s_trend = factor(rep(NA, length(tmb_params$omega_s_trend)))))
  }
  if (is.null(time_varying))
    tmb_map <- c(tmb_map,
      list(b_rw_t = as.factor(matrix(NA, nrow = tmb_data$n_t, ncol = ncol(X_rw_ik)))),
      list(ln_tau_V = as.factor(NA))
    )

  if (nobs_RE[[1]] > 0) tmb_random <- c(tmb_random, "RE")
  if (reml) tmb_random <- c(tmb_random, "b_j")

  if (est_epsilon_model >= 2) {
    # model 2 = re model, model 3 = loglinear-re
    tmb_random <- c(tmb_random, "epsilon_rw")
  }

  if(has_smooths) {
    # random effects for the smooth parameters for p-splines
    tmb_random <- c(tmb_random, "b_smooth")
  }

  if (!is.null(previous_fit)) {
    tmb_params <- previous_fit$tmb_obj$env$parList()
  }

  tmb_data$normalize_in_r <- as.integer(normalize)

  if (!is.null(previous_fit)) tmb_map <- previous_fit$tmb_map
  if (isTRUE(map_rf)) tmb_map <- map_off_rf(tmb_map, tmb_params)
  tmb_map <- c(map, tmb_map)
  if (share_range) tmb_map <- c(tmb_map, list(ln_kappa = factor(c(1L, 1L))))

  for (i in seq_along(start)) {
    message("Initiating ", names(start)[i],
      " at specified starting value of ", start[[i]], ".")
    tmb_params[[names(start)[i]]] <- start[[i]]
  }
  if (length(tmb_params[["ln_kappa"]]) != 2L && "ln_kappa" %in% names(start)) {
    stop("Note that ln_kappa must be a vector of length 2. It should be the same value if `share_range = TRUE`.", call. = FALSE)
  }

  data$sdm_x <- data$sdm_y <- data$sdm_orig_id <- data$sdm_spatial_id <- NULL

  tmb_obj <- TMB::MakeADFun(
    data = tmb_data, parameters = tmb_params, map = tmb_map,
    profile = if (control$profile) "b_j" else NULL,
    random = tmb_random, DLL = "sdmTMB", silent = silent)
  lim <- set_limits(tmb_obj, lower = lower, upper = upper,
    loc = spde$mesh$loc, silent = FALSE)

  out_structure <- structure(list(
    data       = data,
    spde       = spde,
    formula    = original_formula,
    split_formula = split_formula,
    time_varying = time_varying,
    threshold_parameter = thresh$threshold_parameter,
    threshold_function = thresh$threshold_func,
    epsilon_predictor = epsilon_predictor,
    mgcv_mod   = mgcv_mod,
    time       = time,
    family     = family,
    response   = y_i,
    tmb_data   = tmb_data,
    tmb_params = tmb_params,
    tmb_map    = tmb_map,
    tmb_random = tmb_random,
    reml       = reml,
    mgcv       = mgcv,
    lower      = lim$lower,
    upper      = lim$upper,
    priors     = priors,
    nlminb_control = .control,
    control  = control,
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
      h <- stats::optimHess(tmb_opt$par, fn = tmb_obj$fn, gr = tmb_obj$gr,
        upper = lim$upper, lower = lim$lower)
      tmb_opt$par <- tmb_opt$par - solve(h, g)
      tmb_opt$objective <- tmb_obj$fn(tmb_opt$par)
    }
  }
  check_bounds(tmb_opt$par, lim$lower, lim$upper)

  sd_report <- TMB::sdreport(tmb_obj, getJointPrecision = get_joint_precision)
  conv <- get_convergence_diagnostics(sd_report)

  out <- c(out_structure, list(
    model      = tmb_opt,
    tmb_obj    = tmb_obj,
    sd_report  = sd_report,
    gradients  = conv$final_grads,
    bad_eig    = conv$bad_eig))
  `class<-`(out, "sdmTMB")
}

map_off_rf <- function(.map, tmb_params) {
  .map$ln_tau_O <- as.factor(NA)
  .map$ln_tau_E <- as.factor(NA)
  .map$ln_tau_O_trend <- as.factor(NA)
  .map$omega_s_trend <- factor(rep(NA, length(tmb_params$omega_s_trend)))
  .map$ln_kappa <- factor(rep(NA, 2))
  .map$ln_H_input <- factor(rep(NA, 2))
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

# from brms:::rm_wsp()
rm_wsp <- function (x) {
  out <- gsub("[ \t\r\n]+", "", x, perl = TRUE)
  dim(out) <- dim(x)
  out
}
# from brms:::all_terms()
all_terms <- function (x) {
  if (!length(x)) {
    return(character(0))
  }
  if (!inherits(x, "terms")) {
    x <- terms(as.formula(x))
  }
  rm_wsp(attr(x, "term.labels"))
}

get_smooth_terms <- function(terms) {
  x1 <- grep("s\\(", terms)
  x2 <- grep("t2\\(", terms)
  x3 <- grep("by=", terms)
  if (length(x2) > 0) stop("sdmTMB is not set up to work with t2() yet.", call. = FALSE)
  if (length(x3) > 0) warning("s(..., by = ) may not yet be set up correctly within sdmTMB.\n",
    "I.e., it does not match mgcv output exactly.\n",
    "*Use at your own risk.*", call. = FALSE)
  x1
}

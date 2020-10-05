#' @useDynLib sdmTMB, .registration = TRUE
NULL

#' Fit a spatial or spatiotemporal GLMM with TMB
#'
#' Fit a spatial or spatiotemporal GLMM with TMB. Particularly useful for
#' species distribution models and relative abundance index standardization.
#'
#' @param formula Model formula. See the Details section below for how to specify
#'   offsets and threshold parameters. For index standardization, include `0 +
#'   as.factor(year)` (or whatever the time column is called) in the formula.
#' @param data A data frame.
#' @param spde An object from [make_spde()].
#' @param time The time column (as character).
#' @param family The family and link. Supports [gaussian()], [Gamma()],
#'   [binomial()], [poisson()], [sdmTMB::Beta()], [nbinom2()], and [tweedie()].
#' @param time_varying An optional formula describing covariates that should be
#'   modelled as a random walk through time. Leave as `NULL` for a spatial-only
#'   model.
#' @param weights Optional likelihood weights for the conditional model.
#'   Implemented as in \pkg{glmmTMB}. In other words, weights do not have to sum
#'   to one and are not internally modified.
#' @param extra_time Optional extra time slices (e.g., years) to include for
#'   interpolation or forecasting with the predict function. See details section.
#' @param reml Logical: use REML estimation rather than maximum likelihood?
#' @param silent Silent or include optimization details?
#' @param multiphase Logical: estimate the fixed and random effects in phases?
#'   Phases are usually faster and more stable.
#' @param anisotropy Logical: allow for anisotropy? See [plot_anisotropy()].
#' @param control Optimization control options. See [sdmTMBcontrol()].
#' @param enable_priors Should weakly informative priors be enabled?
#'   Experimental and likely for use with the \pkg{tmbstan} package. Note that
#'   the priors are not yet sensible and Jacobian adjustments are not made. If you
#'   are interested in this functionality, please contact the developers.
#' @param ar1_fields Estimate the spatiotemporal random fields as an AR1
#'   process? Note that the parameter `ar1_phi` has been internally bounded
#'   between `-1` and `1` with:  `2 * invlogit(ar1_phi) - 1` i.e. in R `2 *
#'   plogis(ar_phi) - 1`.
#' @param include_spatial Should a separate spatial random field be estimated?
#'   If enabled then there will be separate spatial and spatiotemporal
#'   fields.
#' @param spatial_trend Should a separate spatial field be included in the
#'   trend? Requires spatiotemporal data.
#' @param spatial_only Deprecated; please use `time = NULL` to fit a spatial model.
#' @param nlminb_loops How many times to run [stats::nlminb()] optimization.
#'   Sometimes restarting the optimizer at the previous best values aids
#'   convergence. If the maximum gradient is still too large,
#'   try increasing this to `2`.
#' @param newton_steps How many Newton optimization steps to try with
#'   [stats::optimHess()] after running [stats::nlminb()]. Sometimes aids
#'   convergence.
#' @param mgcv Parse the formula with [mgcv::gam()]?
#' @param previous_fit A previously fitted sdmTMB model to initialize the
#'   optimization with. Can greatly speed up fitting. Note that the data and
#'   model must be set up exactly the same way! However, the `weights` argument
#'   can change, which can be useful for cross-validation.
#' @param quadratic_roots Experimental feature for internal use right now; may
#'   be moved to a branch. Logical: should quadratic roots be calculated? Note:
#'   on the sdmTMB side, the first two coefficients are used to generate the
#'   quadratic parameters. This means that if you want to generate a quadratic
#'   profile for depth, and depth and depth^2 are part of your formula, you need
#'   to make sure these are listed first and that an intercept isn't included.
#'   For example, `formula = cpue ~ 0 + depth + depth2 + as.factor(year)`.
#'
#' @importFrom methods as is
#' @importFrom stats gaussian model.frame model.matrix
#'   model.response terms model.offset
#'
#' @details
#'
#' **Offsets**
#'
#' In the model formula, an offset can be included by including `+ offset` in
#' the model formula (a reserved word). The offset will be included in any
#' prediction. `offset` must be a column in `data`.
#'
#' **Threshold models**
#'
#' A linear break-point relationship for a covariate can be included via `+
#' breakpt(variable)` in the formula, where `variable` is a single covariate
#' corresponding to a column in `data`. In this case the relationship is linear
#' up to a point and then constant.
#'
#' Similarly, a logistic-function threshold model can be included via `+
#' logistic(variable)`. This option models the relationship as a logistic
#' function of the 50% and 95% values. This is similar to length- or size-based
#' selectivity in fisheries, and is parameterized by the points at which f(x) =
#' 0.5 or 0.95. See the vignette.
#'
#' Note that only a single threshold covariate can be included.
#'
#' **Forecasting or interpolating**
#'
#' Extra time slices (e.g., years) can be included for interpolation or
#' forecasting with the predict function via the `extra_time` argument. The
#' predict function requires all time slices to be defined when fitting the
#' model to ensure the various time indices are set up correctly. Be careful if
#' including extra time slices that the model remains identifiable. For example,
#' including `+ as.factor(year)` in `formula` will render a model with no data
#' to inform the expected value in the missing year. [sdmTMB()] makes no attempt
#' to determine if the model makes sense for forecasting or interpolation. The
#' options `time_varying`, `include_spatial`, `ar1_fields`, `time = NULL`
#' provide mechanisms to predict over missing time slices.
#'
#' @export
#'
#' @examples
#' d <- subset(pcod, year >= 2011) # subset for example speed
#' pcod_spde <- make_spde(d, c("X", "Y"), cutoff = 30) # a coarse mesh for example speed
#' plot(pcod_spde)
#'
#' # Tweedie:
#' m <- sdmTMB(density ~ 0 + depth_scaled + depth_scaled2 + as.factor(year),
#'   data = d, time = "year", spde = pcod_spde, family = tweedie(link = "log"))
#' m
#'
#' # Contents of the output object:
#' names(m)
#' m$model
#' r <- m$tmb_obj$report()
#' names(r)
#'
#' # Run extra optimization steps to help convergence:
#' m1 <- run_extra_optimization(m, nlminb_loops = 0, newton_steps = 1)
#' max(m$gradients)
#' max(m1$gradients)
#'
#' # Binomial:
#' pcod_binom <- d
#' pcod_binom$present <- ifelse(pcod_binom$density > 0, 1L, 0L)
#' m_bin <- sdmTMB(present ~ 0 + as.factor(year) + depth_scaled + depth_scaled2,
#'   data = pcod_binom, time = "year", spde = pcod_spde,
#'   family = binomial(link = "logit"))
#'
#' # Gaussian:
#' pcod_gaus <- subset(d, density > 0 & year >= 2013)
#' pcod_spde_gaus <- make_spde(pcod_gaus, c("X", "Y"), cutoff = 30)
#' m_pos <- sdmTMB(log(density) ~ 0 + as.factor(year) + depth_scaled + depth_scaled2,
#'   data = pcod_gaus, time = "year", spde = pcod_spde_gaus)
#'
#' # With splines via mgcv.
#' # Make sure to pre-specify an appropriate basis dimension (`k`) since
#' # the smoothers are not penalized in the current implementation.
#' # See ?mgcv::choose.k
#' m_gam <- sdmTMB(log(density) ~ 0 + as.factor(year) + s(depth_scaled, k = 4),
#'   data = pcod_gaus, time = "year", spde = pcod_spde_gaus)
#'
#' # Fit a spatial only model:
#' m <- sdmTMB(
#'   density ~ depth_scaled + depth_scaled2, data = d,
#'   spde = pcod_spde, family = tweedie(link = "log"))
#'
#' # Spatial-trend example:
#' m <- sdmTMB(density ~ depth_scaled, data = d,
#'   spde = pcod_spde, family = tweedie(link = "log"),
#'   spatial_trend = TRUE, time = "year")
#'
#' r <- m$tmb_obj$report()
#' r$ln_tau_O_trend
#' r$omega_s_trend
#'
#' # Time-varying effects of depth and depth squared:
#' m <- sdmTMB(density ~ 0 + as.factor(year),
#'   time_varying = ~ 0 + depth_scaled + depth_scaled2,
#'   data = d, time = "year", spde = pcod_spde, family = tweedie(link = "log"))
#'
#' # See the b_rw_t estimates; these are the time-varying (random walk) effects.
#' summary(m$sd_report)[1:19,]
#'
#' # Linear breakpoint model on depth:
#' m_pos <- sdmTMB(log(density) ~ 0 + as.factor(year) +
#'     breakpt(depth_scaled) + depth_scaled2, data = pcod_gaus,
#'   time = "year", spde = pcod_spde_gaus)

sdmTMB <- function(formula, data, spde, time = NULL,
  family = gaussian(link = "identity"),
  time_varying = NULL, weights = NULL, extra_time = NULL, reml = FALSE,
  silent = TRUE, multiphase = TRUE, anisotropy = FALSE,
  control = sdmTMBcontrol(), enable_priors = FALSE, ar1_fields = FALSE,
  include_spatial = TRUE, spatial_trend = FALSE,
  spatial_only = "Depreciated",
  nlminb_loops = 1,
  newton_steps = 0,
  mgcv = TRUE,
  previous_fit = NULL,
  quadratic_roots = FALSE) {

  if (!missing(spatial_only)) {
    warning("The argument `spatial_only` is deprecated; please use `time = NULL` instead.",
      call. = FALSE)
    spatial_only <- TRUE
  }
  if (isTRUE(spatial_only) || is.null(time) || identical(length(unique(data[[time]])), 1L)) {
    spatial_only <- TRUE
  } else {
    spatial_only <- FALSE
  }

  if (is.null(time)) {
    time <- "_sdmTMB_time"
    data[[time]] <- 0L
  } else {
    if (sum(is.na(data[[time]])) > 1)
      stop("There is at least one NA value in the time column. ",
        "Please remove it.", call. = FALSE)
  }

  thresh <- check_and_parse_thresh_params(formula, data)
  original_formula <- formula
  formula <- thresh$formula

  if (spatial_trend) {
    numeric_time <- time
    t_i <- as.numeric(data[[numeric_time]])
    t_i <- t_i - mean(unique(t_i), na.rm = TRUE) # middle year = intercept
  } else {
    t_i <- rep(0L, nrow(data))
  }
  contains_offset <- check_offset(formula)

  if (isFALSE(mgcv)) {
    mgcv_mod <- NULL
    X_ij <- model.matrix(formula, data)
    mf <- model.frame(formula, data)
  } else {
    mgcv_mod <- mgcv::gam(formula, data = data) # should be fast enough to not worry
    X_ij <- model.matrix(mgcv_mod)
    mf <- model.frame(mgcv::interpret.gam(formula)$fake.formula, data)
  }

  offset_pos <- grep("^offset$", colnames(X_ij))
  y_i  <- model.response(mf, "numeric")

  if (identical(family$link, "log") && min(y_i, na.rm = TRUE) < 0) {
    stop("`link = 'log'` but the reponse data include values < 0.", call. = FALSE)
  }
  if (identical(family$family, "binomial") && !all(y_i %in% c(0, 1))) {
    stop("`family = 'binomial'` but the reponse data include values other than 0 and 1.", call. = FALSE)
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

  tmb_data <- list(
    y_i        = y_i,
    n_t        = length(unique(data[[time]])),
    t_i        = t_i,
    offset_i   = offset,
    A          = spde$A,
    A_st       = A_st,
    A_spatial_index = data$sdm_spatial_id - 1L,
    year_i     = make_year_i(data[[time]]),
    ar1_fields = as.integer(ar1_fields),
    X_ij       = X_ij,
    X_rw_ik    = X_rw_ik,
    proj_lon   = 0,
    proj_lat   = 0,
    do_predict = 0L,
    calc_se    = 0L,
    pop_pred   = 0L,
    weights_i  = if (!is.null(weights)) weights else rep(1, length(y_i)),
    area_i     = rep(1, length(y_i)),
    normalize_in_r = 0L, # not used
    calc_time_totals = 0L,
    random_walk = !is.null(time_varying),
    enable_priors = as.integer(enable_priors),
    include_spatial = as.integer(include_spatial),
    proj_mesh  = Matrix::Matrix(0, 1, 1), # dummy
    proj_X_ij  = matrix(0, ncol = 1, nrow = 1), # dummy
    proj_X_rw_ik = matrix(0, ncol = 1, nrow = 1), # dummy
    proj_year  = 0, # dummy
    proj_spatial_index = 0, # dummy
    proj_t_i = 0, # dummy
    spde_aniso = make_anisotropy_spde(spde),
    spde       = spde$spde$param.inla[c("M0","M1","M2")],
    anisotropy = as.integer(anisotropy),
    family     = .valid_family[family$family],
    link       = .valid_link[family$link],
    spatial_only = as.integer(spatial_only),
    spatial_trend = as.integer(spatial_trend),
    calc_quadratic_range = as.integer(quadratic_roots),
    X_threshold = thresh$X_threshold,
    proj_X_threshold = 0, # dummy
    threshold_func = thresh$threshold_func
  )
  tmb_data$flag <- 1L # Include data

  b_thresh <- rep(0, 2)
  if (thresh$threshold_func == 2L) b_thresh <- c(0, b_thresh) # logistic

  tmb_params <- list(
    ln_H_input = c(0, 0),
    b_j        = rep(0, ncol(X_ij)),
    ln_tau_O   = 0,
    ln_tau_O_trend = 0,
    ln_tau_E   = 0,
    ln_kappa   = 0,
    thetaf     = 0,
    ln_phi     = 0,
    ln_tau_V   = rep(0, ncol(X_rw_ik)),
    ar1_phi    = 0,
    b_rw_t     = matrix(0, nrow = tmb_data$n_t, ncol = ncol(X_rw_ik)),
    omega_s    = rep(0, n_s),
    omega_s_trend = rep(0, n_s),
    epsilon_st = matrix(0, nrow = n_s, ncol = tmb_data$n_t),
    b_threshold = b_thresh
  )
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
  if (multiphase && is.null(previous_fit)) {
    not_phase1 <- c(tmb_map, list(
      ln_tau_O   = as.factor(NA),
      ln_tau_E   = as.factor(NA),
      ln_tau_V   = factor(rep(NA, ncol(X_rw_ik))),
      ln_tau_O_trend = as.factor(NA),
      omega_s_trend  = factor(rep(NA, length(tmb_params$omega_s_trend))),
      ln_kappa   = as.factor(NA),
      ln_H_input = factor(rep(NA, 2)),
      b_rw_t     = factor(rep(NA, length(tmb_params$b_rw_t))),
      omega_s    = factor(rep(NA, length(tmb_params$omega_s))),
      epsilon_st = factor(rep(NA, length(tmb_params$epsilon_st)))))

    tmb_obj1 <- TMB::MakeADFun(
      data = tmb_data, parameters = tmb_params,
      map = not_phase1, DLL = "sdmTMB", silent = silent)

    tmb_opt1 <- stats::nlminb(
      start = tmb_obj1$par, objective = tmb_obj1$fn,
      gradient = tmb_obj1$gr, control = control)

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
  if (reml) tmb_random <- c(tmb_random, "b_j")

  if (!is.null(previous_fit)) {
    tmb_params <- previous_fit$tmb_obj$env$parList()
  }

  tmb_obj <- TMB::MakeADFun(
    data = tmb_data, parameters = tmb_params, map = tmb_map,
    random = tmb_random, DLL = "sdmTMB", silent = silent)
  # if (tmb_data$normalize_in_r == 1L)
  #   tmb_obj <- TMB::normalize(tmb_obj, flag = "flag")

  if (!is.null(previous_fit)) {
    start <- previous_fit$model$par
  } else {
    start <- tmb_obj$par
  }

  tmb_opt <- stats::nlminb(
    start = start, objective = tmb_obj$fn, gradient = tmb_obj$gr,
    control = control)

  if (nlminb_loops > 1) {
    if(!silent) cat("running extra nlminb loops\n")
    for (i in seq(2, nlminb_loops, length = max(0, nlminb_loops - 1))) {
      temp <- tmb_opt[c("iterations", "evaluations")]
      tmb_opt <- stats::nlminb(
        start = tmb_opt$par, objective = tmb_obj$fn, gradient = tmb_obj$gr,
        control = control)
      tmb_opt[["iterations"]] <- tmb_opt[["iterations"]] + temp[["iterations"]]
      tmb_opt[["evaluations"]] <- tmb_opt[["evaluations"]] + temp[["evaluations"]]
    }
  }
  if (newton_steps > 0) {
    if(!silent) cat("running newtonsteps\n")
    for (i in seq_len(newton_steps)) {
      g <- as.numeric(tmb_obj$gr(tmb_opt$par))
      h <- stats::optimHess(tmb_opt$par, fn = tmb_obj$fn, gr = tmb_obj$gr)
      tmb_opt$par <- tmb_opt$par - solve(h, g)
      tmb_opt$objective <- tmb_obj$fn(tmb_opt$par)
    }
  }

  sd_report <- TMB::sdreport(tmb_obj)
  conv <- get_convergence_diagnostics(sd_report)

  data$sdm_x <- data$sdm_y <- data$sdm_orig_id <- data$sdm_spatial_id <- NULL

  structure(list(
    model      = tmb_opt,
    data       = data,
    spde       = spde,
    formula    = original_formula,
    time_varying = time_varying,
    threshold_parameter = thresh$threshold_parameter,
    threshold_function = thresh$threshold_func,
    mgcv_mod   = mgcv_mod,
    time       = time,
    family     = family,
    response   = y_i,
    tmb_data   = tmb_data,
    tmb_params = tmb_params,
    tmb_map    = tmb_map,
    tmb_random = tmb_random,
    tmb_obj    = tmb_obj,
    reml       = reml,
    mgcv       = mgcv,
    gradients  = conv$final_grads,
    bad_eig    = conv$bad_eig,
    call       = match.call(expand.dots = TRUE),
    sd_report  = sd_report),
    class      = "sdmTMB")
}

set_par_value <- function(opt, par) {
  as.numeric(opt$par[par == names(opt$par)])
}

#' Optimization control options
#'
#' Any arguments to pass to [stats::nlminb()].
#'
#' @param eval.max Maximum number of evaluations of the objective function
#'   allowed.
#' @param iter.max Maximum number of iterations allowed.
#' @param ... Anything else. See the 'Control parameters' section of
#'   [stats::nlminb()].
#'
#' @export
sdmTMBcontrol <- function(eval.max = 1e4, iter.max = 1e4, ...) {
  list(eval.max = eval.max, iter.max = iter.max, ...)
}

get_convergence_diagnostics <- function(sd_report) {
  final_grads <- sd_report$gradient.fixed
  bad_eig <- FALSE
  if (!is.null(sd_report$pdHess)) {
    if (!sd_report$pdHess) {
      warning("The model may not have converged: ",
        "non-positive-definite Hessian matrix.", call. = FALSE)
    } else {
      eigval <- try(1 / eigen(sd_report$cov.fixed)$values, silent = TRUE)
      if (is(eigval, "try-error") || (min(eigval) < .Machine$double.eps * 10)) {
        warning("The model may not have converged: ",
          "extreme or very small eigen values detected.", call. = FALSE)
        bad_eig <- TRUE
      }
      if (any(final_grads > 0.01))
        warning("The model may not have converged. ",
          "Maximum final gradient: ", max(final_grads), ".", call. = FALSE)
    }
  }
  invisible(list(final_grads = final_grads, bad_eig = bad_eig))
}

make_year_i <- function(x) {
  x <- as.integer(as.factor(x))
  x - min(x)
}

check_offset <- function(formula) {
  any(grepl("^offset$",
    gsub(" ", "", unlist(strsplit(as.character(formula), "\\+")))))
}

update_model <- function(object, silent = FALSE) {
  if (!"pop_pred" %in% names(object$tmb_data)) object$tmb_data$pop_pred <- 0L
  if (!"mgcv" %in% names(object)) object$mgcv <- FALSE
  object$tmb_data$weights_i <- rep(1, length(object$tmb_data$y_i))
  object$tmb_data$calc_quadratic_range <- 0L
  object$tmb_data$area_i <- rep(1, length(object$tmb_data$y_i))
  if (!"X_threshold" %in% names(object$tmb_data)) {
    object$tmb_data$X_threshold <- rep(0, nrow(object$data)) # just placeholder
    object$tmb_data$threshold_func <- 0L
    object$tmb_data$proj_X_threshold <- 0 # dummy
    object$tmb_params$b_threshold <- rep(0, 2)
  }
  object$tmb_obj <- TMB::MakeADFun(
    data = object$tmb_data, parameters = object$tmb_params,
    map = object$tmb_map, random = object$tmb_random, DLL = "sdmTMB", silent = silent)
  object
}

check_and_parse_thresh_params <- function(formula, data) {
  terms <- stats::terms(formula)
  terms_labels <- attr(terms, "term.labels")
  if (any(grepl("linear_thresh", terms_labels)) && any(grepl("logistic_thresh", terms_labels))) {
    stop("Please include only a linear (`breakpt`) *or* a logistic threshold.", call. = FALSE)
  }
  if (sum(grepl("linear_thresh", terms_labels)) > 1 || sum(grepl("logistic_thresh", terms_labels)) > 1) {
    stop("Please include only a *single* threshold variable.", call. = FALSE)
  }
  threshold_parameter <- NULL
  if (any(grepl("breakpt", terms_labels))) {
    out <- parse_threshold_formula(formula, "breakpt", terms_labels)
    threshold_parameter <- out$threshold_parameter
    formula <- out$formula
    threshold_function <- "linear"
  }
  if (any(grepl("logistic", terms_labels))) {
    out <- parse_threshold_formula(formula, "logistic", terms_labels)
    threshold_parameter <- out$threshold_parameter
    formula <- out$formula
    threshold_function <- "logistic"
  }
  if(!is.null(threshold_parameter)) {
    if (length(threshold_parameter) > 1) {
      stop("`threshold_parameter` must be a single variable name.", call. = FALSE)
    }
    if (!threshold_parameter %in% names(data)) {
      stop("`threshold_parameter` is not a column in the `data` data frame.", call. = FALSE)
    }
  }

  if (is.null(threshold_parameter)) {
    X_threshold <- rep(0, nrow(data)) # just placeholder
    threshold_func <- 0L
  } else {
    X_threshold <- data[, names(data) == threshold_parameter, drop = TRUE]
    # indexed 1, 2 because 0 will tell TMB not to estimate this:
    threshold_func <- match(threshold_function, c("linear", "logistic"))
  }
  X_threshold <- as.numeric(unlist(X_threshold))
  list(formula = formula, threshold_parameter = threshold_parameter,
    threshold_func = threshold_func, X_threshold = X_threshold
  )
}

parse_threshold_formula <- function(formula, thresh_type_short = "lin_thresh",
  terms_labels) {
  which_thresh <- grep(thresh_type_short, terms_labels)
  temp <- gsub(paste0("^", thresh_type_short, "\\("), "", terms_labels[which_thresh])
  threshold_parameter <- gsub("\\)$", "", temp)
  formula <- stats::update(formula, paste("~ . -", terms_labels[which_thresh]))
  list(formula = formula, threshold_parameter = threshold_parameter)
}


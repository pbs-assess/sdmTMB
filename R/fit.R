#' @useDynLib sdmTMB, .registration = TRUE
NULL

#' Fit a spatial or spatiotemporal GLMM with TMB
#'
#' Fit a spatial or spatiotemporal GLMM with TMB. Particularly useful for
#' species distribution models and relative abundance index standardization.
#'
#' @param formula Model formula. An offset can be included by including `offset`
#'   in the model formula (a reserved word). The offset will be included in any
#'   prediction. For index standardization, include `0 + as.factor(year)` (or
#'   whatever the time column is called) in the formula.
#' @param data A data frame.
#' @param time The time column (as character).
#' @param spde An object from [make_spde()].
#' @param family The family and link. Supports [gaussian()], [Gamma()],
#'   [binomial()], [poisson()], [nbinom2()], and [tweedie()].
#' @param time_varying An optional formula describing covariates that should be
#'   modelled as a random walk through time.
#' @param weights Optional likelihood weights for the conditional model.
#'   Implemented as in \pkg{glmmTMB}. In other words, weights do not have to sum
#'   to one and are not internally modified.
#' @param reml Logical: use REML estimation rather than maximum likelihood?
#' @param silent Silent or include optimization details?
#' @param multiphase Logical: estimate the fixed and random effects in phases?
#'   Phases are usually faster and more stable.
#' @param anisotropy Logical: allow for anisotropy? See [plot_anisotropy()].
#' @param control Optimization control options. See [sdmTMBcontrol()].
#' @param enable_priors Should weakly informative priors be enabled?
#'   Experimental and likely for use with the \pkg{tmbstan} package. Note that
#'   the priors are not yet sensible.
#' @param ar1_fields Estimate the spatiotemporal random fields as an AR1
#'   process? Note that the parameter `ar1_phi` has been internally bounded
#'   between `-1` and `1` with:  `2 * invlogit(ar1_phi) - 1` i.e. in R `2 *
#'   plogis(ar_phi) - 1`.
#' @param include_spatial Should a separate spatial random field be estimated?
#'   If enabled then there will be a separate spatial field and spatiotemporal
#'   fields.
#' @param spatial_trend Should a separate spatial field be included in the
#'   trend? Requires spatiotemporal data.
#' @param normalize Logical: should the normalization of the random effects be
#'   done in R during the outer-optimization step? For some cases, especially
#'   with many knots, this may be faster. In others, it may be slower or suffer
#'   from convergence problems. *Currently disabled!*
#' @param spatial_only Logical: should only a spatial model be fit (i.e. do not
#'   include spatiotemporal random effects)? By default a spatial-only model
#'   will be fit if there is only one unique value in the time column or the
#'   `time` argument is left at its default value of `NULL`.
#' @param quadratic_roots Logical: should quadratic roots be calculated?
#'   Experimental feature for internal use right now. Note: on the sdmTMB side,
#'   the first two coefficients are used to generate the quadratic parameters.
#'   This means that if you want to generate a quadratic profile for depth, and
#'   depth and depth^2 are part of your formula, you need to make sure these are
#'   listed first and that an intercept isn't included. For example, `formula
#'   = cpue ~ 0 + depth + depth2 + as.factor(year)`.
#' @param threshold_parameter Optional parameter to include as a non-linear threshold
#'   relationship. Form can be linear or logistic, and is passed in as a character string,
#'   e.g. "temperature" that is a name of a variable in the data frame
#' @param threshold_function Optional name to include of the threshold function. Defaults to
#'   "linear", in which case a linear breakpoint model is used. Other option is "logistic",
#'   which models the relationship as a function of the 50% and 95% values
#'
#' @importFrom methods as is
#' @importFrom stats gaussian model.frame model.matrix
#'   model.response terms model.offset
#'
#' @export
#'
#' @examples
#' d <- subset(pcod, year >= 2011) # subset for example speed
#' pcod_spde <- make_spde(d$X, d$Y, n_knots = 50) # only 50 knots for example speed
#' plot_spde(pcod_spde)
#'
#' # Tweedie:
#' m <- sdmTMB(density ~ 0 + depth_scaled + depth_scaled2 + as.factor(year),
#' data = d, time = "year", spde = pcod_spde, family = tweedie(link = "log"))
#' m
#'
#' # Contents of the output object:
#' names(m)
#' m$model
#' TMB::sdreport(m$tmb_obj)
#' r <- m$tmb_obj$report()
#' names(r)
#'
#' # Binomial:
#' pcod_binom <- d
#' pcod_binom$present <- ifelse(pcod_binom$density > 0, 1L, 0L)
#' m_bin <- sdmTMB(present ~ 0 + as.factor(year) + depth_scaled + depth_scaled2,
#' data = pcod_binom, time = "year", spde = pcod_spde,
#' family = binomial(link = "logit"))
#'
#' # Gaussian:
#' pcod_gaus <- subset(d, density > 0 & year >= 2013)
#' pcod_spde_gaus <- make_spde(pcod_gaus$X, pcod_gaus$Y, n_knots = 50)
#' m_pos <- sdmTMB(log(density) ~ 0 + as.factor(year) + depth_scaled + depth_scaled2,
#' data = pcod_gaus, time = "year", spde = pcod_spde_gaus)
#'
#' # Fit a spatial only model:
#' m <- sdmTMB(
#' density ~ depth_scaled + depth_scaled2, data = d,
#' spde = pcod_spde, family = tweedie(link = "log"))
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
#' # Experimental calculation of quadratic roots:
#' m <- sdmTMB(density ~ 0 + depth_scaled + depth_scaled2 + as.factor(year),
#'   data = d, time = "year", spde = pcod_spde, family = tweedie(link = "log"),
#'   quadratic_roots = TRUE)
#' .sd_report <- summary(m$sd_report)
#' params <- row.names(.sd_report)
#' params <- .sd_report[grep("quadratic", params), ]
#' params
#' b <- m$model$par[1:2]
#' x <- seq(min(d$depth_scaled), max(d$depth_scaled), length.out = 300)
#' y <- exp(1 + x * b[1] + x^2 * b[2])
#' plot(x, y, type = "l")
#' abline(h = y[which(y == max(y))] * 0.05)
#' abline(v = params[1:2, 1])

sdmTMB <- function(formula, data, time = NULL, spde,
  family = gaussian(link = "identity"),
  time_varying = NULL, weights = NULL, reml = FALSE,
  silent = TRUE, multiphase = TRUE, anisotropy = FALSE,
  control = sdmTMBcontrol(), enable_priors = FALSE, ar1_fields = FALSE,
  include_spatial = TRUE, spatial_trend = FALSE,
  normalize = FALSE,
  spatial_only = identical(length(unique(data[[time]])), 1L),
  quadratic_roots = FALSE, threshold_parameter = NULL, threshold_function="linear") {

  if (isTRUE(normalize)) {
    warning("`normalize` is currently disabled and doesn't do anything.")
    normalize <- FALSE
  }

  if (is.null(time)) {
    time <- "_sdmTMB_time"
    data[[time]] <- 0L
  } else {
    if (sum(is.na(data[[time]])) > 1)
      stop("There is at least one NA value in the time column. ",
        "Please remove it.", call. = FALSE)
  }

  if(!is.null(threshold_parameter)) {
    if(threshold_function %in% c("linear","logistic") == FALSE) {
      stop("If you're trying to specify a threshold model, please specify the function type as 'linear' or 'logistic'.", call. = FALSE)
    }
    if(threshold_parameter %in% names(data) == FALSE) {
      stop("If you're trying to specify a threshold model, make sure the name of the variable is in the data frame.", call. = FALSE)
    }
  }

  if (spatial_trend) {
    numeric_time <- time
    t_i <- as.numeric(data[[numeric_time]])
    t_i <- t_i - mean(unique(t_i), na.rm = TRUE) # middle year = intercept
  } else {
    t_i <- rep(0L, nrow(data))
  }
  contains_offset <- check_offset(formula)
  # X_ij contains linear fixed effects. Don't include non-linear threshold parameter
  # in it if it's passed in
  if(is.null(threshold_parameter)) {
    X_ij <- model.matrix(formula, data)
    mf   <- model.frame(formula, data)
    X_threshold = rep(0, nrow(X_ij)) # just placeholder
    threshold_func = 0
  } else {
    # omit threshold_parameter from linear fixed effects model matrix / frame
    X_ij <- model.matrix(formula, data[,-which(names(data) == threshold_parameter)])
    mf   <- model.frame(formula, data[,-which(names(data) == threshold_parameter)])
    X_threshold = data[,which(names(data) == threshold_parameter)]
    # indexed 1-2 because 0 will tell TMB not to estimate this
    threshold_func = as.numeric(match(threshold_function,c("linear","logistic")))
  }
  offset_pos <- grep("^offset$", colnames(X_ij))
  y_i  <- model.response(mf, "numeric")
  offset <- as.vector(model.offset(mf))
  if (is.null(offset)) offset <- rep(0, length(y_i))

  if (!is.null(time_varying)) {
    X_rw_ik <- model.matrix(time_varying, data)
  } else {
    X_rw_ik <- matrix(0, nrow = nrow(data), ncol = 1)
  }
  # Stuff needed for spatiotemporal A matrix:
  data$sdm_orig_id <- seq(1, nrow(data))
  data$sdm_x <- spde$x
  data$sdm_y <- spde$y
  fake_data <- unique(data.frame(sdm_x = spde$x, sdm_y = spde$y))
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
    weights_i  = if (!is.null(weights)) weights else rep(1, length(y_i)),
    area_i     = rep(1, length(y_i)),
    normalize_in_r = as.integer(normalize),
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
    X_threshold = as.numeric(unlist(X_threshold)),
    proj_X_threshold = 0, # dummy
    threshold_func = as.integer(threshold_func)
  )
  tmb_data$flag <- 1L # Include data

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
    b_threshold = rep(0,4)
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

  if(is.null(threshold_parameter)) {
    tmb_map <- c(tmb_map, list(b_threshold = factor(rep(NA, 3))))
  }
  if (multiphase) {
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

  tmb_obj <- TMB::MakeADFun(
    data = tmb_data, parameters = tmb_params, map = tmb_map,
    random = tmb_random, DLL = "sdmTMB", silent = silent)
  if (tmb_data$normalize_in_r == 1L)
    tmb_obj <- TMB::normalize(tmb_obj, flag = "flag")

  tmb_opt <- stats::nlminb(
    start = tmb_obj$par, objective = tmb_obj$fn, gradient = tmb_obj$gr,
    control = control)

  sd_report <- TMB::sdreport(tmb_obj)
  conv <- get_convergence_diagnostics(sd_report)

  data$sdm_x <- data$sdm_y <- data$sdm_orig_id <- data$sdm_spatial_id <- NULL

  structure(list(
    model      = tmb_opt,
    data       = data,
    spde       = spde,
    formula    = formula,
    time_varying = time_varying,
    threshold_parameter = threshold_parameter,
    threshold_function = threshold_function,
    time       = time,
    family     = family,
    response   = y_i,
    tmb_data   = tmb_data,
    tmb_params = tmb_params,
    tmb_map    = tmb_map,
    tmb_random = tmb_random,
    tmb_obj    = tmb_obj,
    reml       = reml,
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
  object$tmb_data$weights_i <- rep(1, length(object$tmb_data$y_i))
  object$tmb_data$calc_quadratic_range <- 0L
  object$tmb_data$area_i <- rep(1, length(object$tmb_data$y_i))
  object$tmb_obj <- TMB::MakeADFun(
    data = object$tmb_data, parameters = object$tmb_params,
    map = object$tmb_map, random = object$tmb_random, DLL = "sdmTMB", silent = silent)
  object
}


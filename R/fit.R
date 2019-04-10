#' @useDynLib sdmTMB
NULL

#' Fit a spatiotemporal GLMM with TMB, e.g. for a species distribution model.
#'
#' @param formula Model formula. For index standardization you will want to
#'   include `0 + as.factor(your_time_column)`.
#' @param data A data frame.
#' @param time The time column (as character).
#' @param spde An object from [make_spde()].
#' @param family The family and link. Supports [gaussian()], [Gamma()],
#'   [binomial()], [poisson()], [nbinom2()], and [tweedie()].
#' @param time_varying An optional formula describing covariates that should be
#'   modelled as a random walk through time.
#' @param silent Silent or optimization details?
#' @param multiphase Estimate the fixed and random effects in phases for speed?
#' @param anisotropy Logical: allow for anisotropy?
#' @param control Optimization control options. See [sdmTMBcontrol()].
#' @param enable_priors Should weakly informative priors be enabled?
#'   (experimental and likely for use with the \pkg{tmbstan} package)
#' @param ar1_fields Estimate the spatiotemporal random fields as an AR1
#'   process? Note that the parameter `ar1_phi` has been internally bounded
#'   between -1 and 1 with:  `2 * invlogit(ar1_phi) - 1` i.e. in R ` 2 *
#'   plogis(ar_phi) - 1`.
#' @param include_spatial Should a separate spatial random field the estimated?
#'   If enabled then there will be a separate spatial field and spatiotemporal
#'   fields.
#' @param spatial_trend Should a separate spatial field be included in the
#'   trend? This works if hauls can be viewed as replicates of grid cell
#'   observations, and only when other spatiotemporal components are not
#'   estimated.
#' @param normalize Logical: should the normalization of the random effects
#'   be done in R during the outer-optimization step? For some cases,
#'   especially with many knots, this may be faster. In others, it may be slower
#'   or suffer from convergence problems.
#' @param spatial_only Logical: should only a spatial model be fit (i.e. do not
#'   include spatiotemporal random effects)? By default a spatial-only model
#'   will be fit if there is only one unique value in the time column or the
#'   `time` argument is left at its default value of `NULL`.
#'
#' @importFrom methods as is
#' @importFrom stats gaussian model.frame model.matrix
#'   model.response terms
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

sdmTMB <- function(data, formula, time = NULL, spde, family = gaussian(link = "identity"),
  time_varying = NULL, silent = TRUE, multiphase = TRUE, anisotropy = FALSE,
  control = sdmTMBcontrol(), enable_priors = FALSE, ar1_fields = FALSE,
  include_spatial = TRUE, spatial_trend = FALSE, normalize = FALSE,
  spatial_only = identical(length(unique(data[[time]])), 1L)) {

  if (is.null(time)) {
    time <- "_sdmTMB_time"
    data[[time]] <- 0L
  } else {
    if (sum(is.na(data[[time]])) > 1)
      stop("There is at least one NA value in the time column. ",
        "Please remove it.", call. = FALSE)
  }

  if (spatial_trend) {
    numeric_time <- time
    t_i <- as.numeric(data[[numeric_time]])
    t_i <- t_i - mean(unique(t_i), na.rm = TRUE) # middle year = intercept
  } else {
    t_i <- rep(0L, nrow(data))
  }

  X_ij <- model.matrix(formula, data)
  mf   <- model.frame(formula, data)
  y_i  <- model.response(mf, "numeric")

  if (!is.null(time_varying))
    X_rw_ik <- model.matrix(time_varying, data)
  else
    X_rw_ik <- matrix(0, nrow = nrow(data), ncol = 1)

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
    spatial_trend = as.integer(spatial_trend)
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
    omega_s_trend    = rep(0, n_s),
    epsilon_st = matrix(0, nrow = n_s, ncol = tmb_data$n_t)
  )

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
    tmb_params$b_j <- set_par_value(tmb_opt1, "b_j")
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

  args <- get_args()
  structure(list(
    model      = tmb_opt,
    data       = data,
    spde       = spde,
    formula    = formula,
    time_varying = time_varying,
    time       = time,
    family     = family,
    response   = y_i,
    tmb_data   = tmb_data,
    tmb_params = tmb_params,
    tmb_map    = tmb_map,
    tmb_random = tmb_random,
    tmb_obj    = tmb_obj,
    gradients  = conv$final_grads,
    bad_eig    = conv$bad_eig,
    sd_report  = sd_report,
    args       = args),
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

# Get a list of the arguments used within any function call
get_args <- function(){
  def.call <- sys.call(-1)
  def <- get(as.character(def.call[[1]]), mode="function", sys.frame(-2))
  act.call <- match.call(definition = def, call = def.call)
  def <- as.list(def)
  def <- def[-length(def)]
  act <- as.list(act.call)[-1]

  def.nm <- names(def)
  act.nm <- names(act)
  inds <- def.nm %in% act.nm
  out <- def
  out[inds] <- act
  out
}

make_year_i <- function(x) {
  x <- as.integer(as.factor(x))
  x - min(x)
}

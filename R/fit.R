# library(INLA)
# set.seed(1)
# n = 20
# x_loc = cbind(runif(n), runif(n))
# mesh1 = inla.mesh.2d(loc = x_loc, max.edge = 0.1)
# mesh2 = inla.mesh.2d(loc = x_loc, max.edge = c(0.1, 0.4))
# mesh3 = inla.mesh.2d(loc.domain = x_loc, max.edge = c(0.1, 0.4), )
# plot(mesh3)
# points(x_loc, col = "red", pch = 20, cex = 2)
#
# loc_xy <- data.frame(pcod$X, pcod$Y)
# bnd = inla.nonconvex.hull(as.matrix(loc_xy), convex = -0.05)
# mesh = inla.mesh.2d(
#   boundary = bnd,
#   max.edge = c(20, 50),
#   offset = -0.05,
#   cutoff = c(2, 5),
#   min.angle = 10
# )
# plot(mesh)
# points(loc_xy, col = "red", pch = 20, cex = 1)
#
# n_knots <- 300
# knots <- stats::kmeans(x = loc_xy, centers = n_knots)
# loc_centers <- knots$centers
# # loc_xy <- cbind(loc_xy, cluster = knots$cluster, time = as.numeric(time))
# mesh <- INLA::inla.mesh.create(loc_centers, refine = TRUE)
# bnd = inla.nonconvex.hull(as.matrix(loc_centers), convex = -0.10)
#
# mesh = inla.mesh.2d(
#   boundary = bnd,
#   max.edge = c(5, 20),
#   offset = c(-0.1, -0.1),
#   cutoff = c(10, 20),
#   min.angle = c(10, 10)
# )
#
# plot(mesh, main = NA, edge.color = "grey60")
# # points(pcod$X, pcod$Y, pch = 21, col = "#00000070")
# points(loc_centers, pch = 20, col = "red")


#' @useDynLib sdmTMB
NULL

#' Construct an SPDE mesh
#'
#' @param x X numeric vector.
#' @param y Y numeric vector.
#' @param n_knots The number of knots.
#' @param seed Random seed. Affects [stats::kmeans()] determination of knot locations.
#'
#' @importFrom graphics points
#' @export
#' @examples
#' sp <- make_spde(pcod$X, pcod$Y, n_knots = 25)
#' plot_spde(sp)
make_spde <- function(x, y, n_knots, seed = 42) {
  loc_xy <- cbind(x, y)
  if (n_knots >= nrow(loc_xy)) {
    warning("Reducing `n_knots` to be one less than the ",
      "number of data points.")
    n_knots <- nrow(loc_xy) - 1
  }
  set.seed(seed)
  knots <- stats::kmeans(x = loc_xy, centers = n_knots)
  loc_centers <- knots$centers
  mesh <- INLA::inla.mesh.create(loc_centers, refine = TRUE)
  spde <- INLA::inla.spde2.matern(mesh)
  list(x = x, y = y, mesh = mesh, spde = spde, cluster = knots$cluster, loc_centers = loc_centers)
}

#' @param object Output from [make_spde()].
#' @rdname make_spde
#' @importFrom graphics plot
#' @export
plot_spde <- function(object) {
  plot(object$mesh, main = NA, edge.color = "grey60", asp = 1)
  points(object$x, object$y, pch = 21, col = "#00000070")
  points(object$loc_centers, pch = 20, col = "red")
}

# from TMB examples repository:
make_anisotropy_spde <- function(spde) {
  inla_mesh <- spde$mesh
  Dset <- 1:2
  inla_mesh <- spde$mesh
  # Triangle info
  TV <- inla_mesh$graph$tv # Triangle to vertex indexing
  V0 <- inla_mesh$loc[TV[, 1], Dset] # V = vertices for each triangle
  V1 <- inla_mesh$loc[TV[, 2], Dset]
  V2 <- inla_mesh$loc[TV[, 3], Dset]
  E0 <- V2 - V1 # E = edge for each triangle
  E1 <- V0 - V2
  E2 <- V1 - V0
  # Calculate Areas
  TmpFn <- function(Vec1, Vec2) abs(det(rbind(Vec1, Vec2)))
  Tri_Area <- rep(NA, nrow(E0))
  for (i in 1:length(Tri_Area)) Tri_Area[i] <- TmpFn(E0[i, ], E1[i, ]) / 2

  list(
    n_s      = spde$spde$n.spde,
    n_tri    = nrow(TV),
    Tri_Area = Tri_Area,
    E0       = E0,
    E1       = E1,
    E2       = E2,
    TV       = TV - 1,
    G0       = spde$spde$param.inla$M0,
    G0_inv   = as(Matrix::diag(1/Matrix::diag(spde$spde$param.inla$M0)), "dgTMatrix"))
}

#' Fit a spatiotemporal GLMM with TMB, e.g. for a species distribution model.
#'
#' @param data A data frame.
#' @param formula Model formula. For index standardization you will want to
#'   include `0 + as.factor(your_time_column)`.
#' @param time The time column (as character).
#' @param spde An object from [make_spde()].
#' @param family The family and link. Supports [gaussian()], [Gamma()],
#'   [binomial()], [poisson()], [nbinom2()], and [tweedie()].
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
#'
#' @importFrom methods as
#' @importFrom stats gaussian model.frame model.matrix
#'   model.response terms
#'
#' @export
#'
#' @examples
#' d <- subset(pcod, year >= 2011) # subset for example speed
#' pcod_spde <- make_spde(d$X, d$Y, n_knots = 100) # only 100 knots for example speed
#' plot_spde(pcod_spde)
#'
#' # Tweedie:
#' m <- sdmTMB(
#'   d, density ~ 0 + as.factor(year) + depth_scaled + depth_scaled2,
#'   time = "year", spde = pcod_spde, family = tweedie(link = "log"),
#'   silent = FALSE)
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
#' m_bin <- sdmTMB(pcod_binom,
#'   present ~ 0 + as.factor(year) + depth_scaled + depth_scaled2,
#'   time = "year", spde = pcod_spde, family = binomial(link = "logit"))
#'
#' # Gaussian:
#' pcod_gaus <- subset(d, density > 0 & year >= 2013)
#' pcod_spde_gaus <- make_spde(pcod_gaus$X, pcod_gaus$Y, n_knots = 50)
#' m_pos <- sdmTMB(pcod_gaus,
#'   log(density) ~ 0 + as.factor(year) + depth_scaled + depth_scaled2,
#'   time = "year", spde = pcod_spde_gaus)
#'
#' \dontrun{
#' # Stan sampling (warning: slow going and priors are flat).
#'
#' # Must load tmbstan first and then TMB and/or sdmTMB
#' # or you will get the error `"is_Null_NS" not resolved from current
#' # namespace (rstan)`
#' # Restart R session, then:
#' library(tmbstan)
#'
#' # Then:
#' library(sdmTMB)
#'
#' # Then:
#' set.seed(42)
#' pcod_pos <- subset(pcod, year > 2013 & density > 0)
#' pcod_pos_spde <- make_spde(pcod_pos$X/10, pcod_pos$Y/10, n_knots = 200) # scale UTMs for #' Stan
#' m <- sdmTMB(pcod_pos,
#'  log(density) ~ 0 + as.factor(year) + depth_scaled + depth_scaled2,
#'  time = "year", spde = pcod_pos_spde)
#' m_stan <- tmbstan(m$tmb_obj, chains = 1, iter = 200, cores=1,
#'   init = "last.par.best", control = list(adapt_delta = 0.80, max_treedepth = 20),
#'   seed = 123, laplace = T)
#'
#' pars <- c('b_j', 'ln_tau_O', 'ln_tau_E', 'ln_kappa', 'ln_phi')
#' m_stan2 <- tmbstan(m$tmb_obj, chains = 1, iter = 200, cores=1,
#'   init = "last.par.best", control = list(adapt_delta = 0.80, max_treedepth = 20),
#'   seed = 123, laplace = F, pars = pars)
#'
#' m_stan
#' }

sdmTMB <- function(data, formula, time, spde, family = gaussian(link = "identity"),
  time_varying = NULL,
  silent = TRUE, multiphase = TRUE, anisotropy = FALSE, control = sdmTMBcontrol(),
  enable_priors = FALSE, ar1_fields = FALSE, include_spatial = TRUE) {

  X_ij <- model.matrix(formula, data)
  mf   <- model.frame(formula, data)
  y_i  <- model.response(mf, "numeric")

  if (!is.null(time_varying))
    X_rw_ik <- model.matrix(time_varying, data)
  else
    X_rw_ik <- matrix(0, nrow = nrow(data), ncol = 1)

  spatial_only <- identical(length(unique(data[[time]])), 1L)

  tmb_data <- list(
    y_i        = y_i,
    n_t        = length(unique(data[[time]])),
    n_s        = nrow(spde$mesh$loc),
    s_i        = spde$cluster - 1L,
    year_i     = as.numeric(as.factor(as.character(data[[time]]))) - 1L,
    year_prev_i= as.numeric(as.factor(as.character(data[[time]]))) - 2L,
    ar1_fields = as.integer(ar1_fields),
    X_ij       = X_ij,
    X_rw_ik    = X_rw_ik,
    proj_lon    = 0,
    proj_lat    = 0,
    do_predict = 0L,
    calc_se    = 0L,
    calc_time_totals = 0L,
    random_walk = !is.null(time_varying),
    enable_priors = as.integer(enable_priors),
    include_spatial = as.integer(include_spatial),
    proj_mesh  = Matrix::Matrix(0, 1, 1), # dummy
    proj_X_ij  = matrix(0, ncol = 1, nrow = 1), # dummy
    proj_X_rw_ik = matrix(0, ncol = 1, nrow = 1), # dummy
    proj_year  = 0, # dummy
    spde_aniso = make_anisotropy_spde(spde),
    spde       = spde$spde$param.inla[c("M0","M1","M2")],
    anisotropy = as.integer(anisotropy),
    family     = .valid_family[family$family],
    link       = .valid_link[family$link],
    spatial_only = as.integer(spatial_only)
  )

  tmb_params <- list(
    ln_H_input = c(0, 0),
    b_j        = rep(0, ncol(X_ij)),
    ln_tau_O   = 0,
    ln_tau_E   = 0,
    ln_kappa   = 0,
    thetaf     = 0,
    ln_phi     = 0,
    ln_tau_V   = rep(0, ncol(X_rw_ik)),
    ar1_phi    = 0,
    b_rw_t     = matrix(0, nrow = tmb_data$n_t, ncol = ncol(X_rw_ik)),
    omega_s    = rep(0, tmb_data$n_s),
    epsilon_st = matrix(0, nrow = tmb_data$n_s, ncol = tmb_data$n_t)
  )

  # Mapping off params as needed:
  tmb_map <- list()
  if (!anisotropy)
    tmb_map <- c(tmb_map, list(ln_H_input = factor(rep(NA, 2))))
  if (!ar1_fields)
    tmb_map <- c(tmb_map, list(ar1_phi = as.factor(NA)))
  if (family$family == "binomial")
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
  if (!is.null(time_varying)) tmb_random <- c(tmb_random, "b_rw_t")

  if (!include_spatial) {
    tmb_map <- c(tmb_map, list(
      ln_tau_O = as.factor(NA),
      omega_s  = factor(rep(NA, length(tmb_params$omega_s)))))
  }

  if (is.null(time_varying))
    tmb_map <- c(tmb_map,
      list(b_rw_t = as.factor(matrix(NA, nrow = tmb_data$n_t, ncol = ncol(X_rw_ik)))),
      list(ln_tau_V = as.factor(NA))
    )

  tmb_obj <- TMB::MakeADFun(
    data = tmb_data, parameters = tmb_params, map = tmb_map,
    random = tmb_random, DLL = "sdmTMB", silent = silent)

  tmb_opt <- stats::nlminb(
    start = tmb_obj$par, objective = tmb_obj$fn, gradient = tmb_obj$gr,
    control = control)

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
      tmb_obj    = tmb_obj),
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

#' Tweedie family
#'
#' @param link The link.
#' @export
#' @examples
#' tweedie(link = "log")
tweedie <- function(link = "log") {
  linktemp <- substitute(link)
  if (!is.character(linktemp))
    linktemp <- deparse(linktemp)
  okLinks <- c("inverse", "log", "identity")
  if (linktemp %in% okLinks)
    stats <- stats::make.link(linktemp)
  else if (is.character(link))
    stats <- stats::make.link(link)

  list(family = "tweedie", link = linktemp, linkfun = stats$linkfun,
    linkinv = stats$linkinv)
}

#' Negative binomial family
#'
#' Specifically the NB2 parameterization where the variance grows quadratically
#' with the mean.
#'
#' @param link The link. Must be `"log"`.
#' @export
#' @examples
#' nbinom2(link = "log")
nbinom2 <- function(link = "log") {
  if (!identical("log", as.character(substitute(link))))
    stop("Link must be 'log' for this implementation of the nbinom2 family.",
      call. = FALSE)

  linktemp <- substitute(link)
  if (!is.character(linktemp))
    linktemp <- deparse(linktemp)
  okLinks <- c("log")
  if (linktemp %in% okLinks)
    stats <- stats::make.link(linktemp)
  else if (is.character(link))
    stats <- stats::make.link(link)

  list(family = "nbinom2", link = linktemp, linkfun = stats$linkfun,
    linkinv = stats$linkinv)
}


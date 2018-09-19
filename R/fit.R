#' @useDynLib sdmTMB
#' @importFrom TMB MakeADFun
NULL

#' Construct an SPDE mesh
#'
#' @param x X numeric vector.
#' @param y Y numeric vector.
#' @param n_knots The number of knots.
#' @param plot Logical: should a plot be made?
#'
#' @importFrom graphics points
#' @export
#' @examples
#' make_spde(pcod$X, pcod$Y, n_knots = 25)
make_spde <- function(x, y, n_knots, plot = FALSE) {
  loc_xy <- data.frame(x, y)
  knots <- stats::kmeans(x = loc_xy, centers = n_knots)
  loc_centers <- knots$centers
  # loc_xy <- cbind(loc_xy, cluster = knots$cluster, time = as.numeric(time))
  mesh <- INLA::inla.mesh.create(loc_centers, refine = TRUE)
  spde <- INLA::inla.spde2.matern(mesh)
  if (plot) {
    plot(mesh, main = NA, edge.color = "grey60")
    points(x, y, pch = 21, col = "#00000070")
    points(loc_centers, pch = 20, col = "red")
  }
  list(mesh = mesh, spde = spde, cluster = knots$cluster, loc_centers = loc_centers)
}

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
    G0_inv   = as(Matrix::diag(1/Matrix::diag(spde$spde$param.inla$M0)),
      "dgTMatrix"))
}

#' Tweedie family
#'
#' @param link The link (must be log)
#' @export
#' @importFrom assertthat assert_that
#' @examples
#' tweedie(link = "log")
tweedie <- function(link = "log") {
  assert_that(identical("log", as.character(substitute(link))),
    msg = "Link must be 'log' for the tweedie family.")
  list(family = "tweedie", link = "log")
}

# Check a family object
#
# @param family The family
check_family <- function(family) {
  list(family = family$family, link = family$link)
}


#' sdmTMB
#'
#' @param data A data frame.
#' @param formula Model formula. For index standardization you will want to
#'   include `0 + as.factor(your_time_column)`.
#' @param time The time column (as character).
#' @param spde An object from [make_spde()].
#' @param family The family and link.
#' @param silent Silent or optimization details?
#' @param multiphase Estimate the fixed and random effects in phases for speed?
#'
#' @importFrom methods as
#' @importFrom stats gaussian model.frame model.matrix
#'   model.response terms
#'
#' @export
#'
#' @examples
#' pcod_spde <- make_spde(pcod$X, pcod$Y, n_knots = 30)
#' m <- sdmTMB(
#'   pcod, density ~ 0 + as.factor(year) + depth_scaled + depth_scaled2,
#'   time = "year", spde = pcod_spde, family = tweedie(link = "log")
#' )

sdmTMB <- function(data, formula, time, spde, family = gaussian(link = "log"),
  silent = FALSE, multiphase = TRUE) {

  X_ij <- model.matrix(formula, data)
  mf   <- model.frame(formula, data)
  y_i  <- model.response(mf, "numeric")

  proj_mesh <- m <- Matrix::Matrix(0, 1, 1) # dummy
  proj_X_ij <- matrix(0, ncol = 1, nrow = 1) # dummy

  family <- check_family(family)
  family_integer <- switch(family$family,
    gaussian = 1L,
    tweedie  = 2L,
    binomial = 3L,
    stop("Family not implemented.")
  )
  link_integer <- switch(family$link,
    identity = 1L,
    log      = 2L,
    logit    = 3L,
    stop("Link not implemented.")
  )

  tmb_data <- list(
    y_i = y_i,
    n_t = length(unique(data[[time]])),
    n_s = nrow(spde$mesh$loc),
    s_i = spde$cluster - 1L,
    year_i = as.numeric(as.factor(as.character(data[[time]]))) - 1L,
    X_ij = X_ij,
    do_predict = 0L,
    proj_mesh = proj_mesh,
    proj_X_ij = proj_X_ij,
    spde = make_anisotropy_spde(spde),
    family = family_integer,
    link = link_integer
  )

  tmb_params <- list(
    ln_H_input = c(-0.5, -0.5),
    b_j        = rep(0, ncol(X_ij)),
    ln_tau_E   = 0.5,
    ln_kappa   = -1,
    logit_p    = 0,
    log_phi    = 2,
    epsilon_st = matrix(0, nrow = tmb_data$n_s, ncol = tmb_data$n_t)
  )

  if (multiphase) {
    not_phase1 <- list(
      ln_tau_E   = as.factor(NA),
      ln_kappa   = as.factor(NA),
      ln_H_input = factor(rep(NA, 2)),
      epsilon_st = factor(rep(NA, length(tmb_params$epsilon_st))) )

    tmb_obj1 <- TMB::MakeADFun(
      data = tmb_data, parameters = tmb_params, map = not_phase1, DLL = "sdmTMB",
      silent = silent)

    tmb_opt1 <- stats::nlminb(
      start = tmb_obj1$par, objective = tmb_obj1$fn, gradient = tmb_obj1$gr,
      control = list(eval.max = 1e4, iter.max = 1e4)
    )

    tmb_params$b_j     <- as.numeric(tmb_opt1$par["b_j"     == names(tmb_opt1$par)])
    tmb_params$logit_p <- as.numeric(tmb_opt1$par["logit_p" == names(tmb_opt1$par)])
    tmb_params$log_phi <- as.numeric(tmb_opt1$par["log_phi" == names(tmb_opt1$par)])
  }

  tmb_random <- c("epsilon_st")
  tmb_obj <- TMB::MakeADFun(
    data = tmb_data, parameters = tmb_params,
    random = tmb_random, DLL = "sdmTMB", silent = silent
  )

  tmb_opt <- stats::nlminb(
    start = tmb_obj$par, objective = tmb_obj$fn, gradient = tmb_obj$gr,
    control = list(eval.max = 1e4, iter.max = 1e4)
  )

  structure(list(
    model = tmb_opt,
    data = data,
    spde = spde,
    formula = formula,
    time = time,
    family = family,
    tmb_data = tmb_data,
    tmb_params = tmb_params,
    response = y_i,
    random = tmb_random,
    tmb_obj = tmb_obj),
    class = "sdmTMB")
}

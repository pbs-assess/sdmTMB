set.seed(1)
predictor_dat <- data.frame(
  X = runif(300), Y = runif(300),
  b1 = rnorm(300), b2 = rnorm(300),
  year = rep(1:6, each = 50)
)
mesh <- make_mesh(predictor_dat, xy_cols = c("X", "Y"), cutoff = 0.1)

sim_dat <- sdmTMB_simulate(
  formula = ~ 1 + b1 + b2,
  data = predictor_dat,
  time = "year",
  mesh = mesh,
  family = gaussian(),
  range = 0.5,
  sigma_E = 0.2,
  sigma_O = 0.7,
  phi = 0.4,
  seed = 123,
  B = c(0.2, -0.4, 0.3)
)
fit <- sdmTMB(observed ~ 1 + b1 + b2, data = sim_dat, mesh = mesh, time = "year")

fixef <- function(x) {
  b <- tidy(fit)
  stats::setNames(b$estimate, b$term)
}

fixef(fit)
fe <- fixef(fit)
X <- fit$tmb_data$X_ij
VarF <- var(as.vector(fixef(fit) %*% t(X))) # variance from fixed-effects
b <- tidy(fit, "ran_par")
sigma_O <- b$estimate[b$term == "sigma_O"] # spatial standard deviation
sigma_E <- b$estimate[b$term == "sigma_E"] # spatiotemporal standard deviation
VarO <- sigma_O^2
VarE <- sigma_E^2
VarR <- var(as.vector(residuals(fit))) # residual variance
VarF/(VarF + VarO + VarE + VarR)
VarO/(VarF + VarO + VarE + VarR)
VarE/(VarF + VarO + VarE + VarR)

# https://github.com/glmmTMB/glmmTMB/blob/master/glmmTMB/inst/misc/rsqglmm.R

r2.sdmTMB <- function(x, which_fixef = NULL) {

  if (!is(fit, "sdmTMB")) {
    stop("x must be a model of class sdmTMB.", call. = FALSE)
  }
  if (!identical(x$family$family, "gaussian"))
    stop("r2.sdmTMB() currently only works for Gaussian models.", call. = FALSE)
  if (!is.null(x$spatial_varying)) {
    stop("r2.sdmTMB() currently does not work with spatially varying coefficient models.", call. = FALSE)
  }

  fe <- fixef(fit)
  X <- fit$tmb_data$X_ij
  if (!is.null(which_fixef)) {
    assert_that(max(which_fixef) <= length(fe))
    assert_that(min(which_fixef) >= 1)
    assert_that(is.numeric(which_fixef))
    assert_that(all(which_fixef %in% seq_along(fe)))
    message("Including fixed effects: ", paste(names(fe)[which_fixef], collapse = ", "))
    fe <- fe[which_fixef]
    X <- X[,which_fixef,drop=FALSE]
  }
  VarF <- var(as.vector(fe %*% t(X))) # variance from fixed-effects

  b <- tidy(fit, "ran_par")

  if (x$tmb_data$include_spatial == 1L) {
    varO <- b$estimate[b$term == "sigma_O"]^2 # spatial variance
  } else {
    varO <- 0
  }

  if (x$tmb_data$spatial_only == 0L) {
    varE <- b$estimate[b$term == "sigma_E"]^2 # spatiotemporal variance
  } else {
    varE <- 0
  }

  VarR <- suppressMessages(var(as.vector(residuals(fit)))) # residual variance

  marg_r2 <- VarF/(VarF + VarO + VarE + VarR)

  if (varO != 0) {
    cond_rf_sp_r2 <- (VarO + VarF)/(VarF + VarO + VarE + VarR)
  } else {
    cond_rf_sp_r2 <- NULL
  }
  if (varE != 0) {
    cond_rf_spt_r2 <- (VarE + VarF)/(VarF + VarO + VarE + VarR)
  } else {
    cond_rf_sp_r2 <- NULL
  }

  if (varE != 0 || varO != 0) {
    cond_rf_r2 <- (VarE + VarO + VarF)/(VarF + VarO + VarE + VarR)
  } else {
    cond_rf_r2 <- NULL
  }

  list(
    marg = marg_r2,
    cond_rf_sp = cond_rf_sp_r2,
    cond_rf_spt = cond_rf_spt_r2,
    cond_rf_all = cond_rf_r2
  )
}

r2.sdmTMB(fit)

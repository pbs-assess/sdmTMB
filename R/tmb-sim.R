#' Simulate from a spatial/spatiotemporal model
#'
#' @description
#' `r lifecycle::badge("experimental")`
#'
#' This version (vs. [sdmTMB_sim()]) (1) uses TMB for simulation and is
#' therefore **much** faster and more flexible, (2) is set up to take a formula and
#' a data frame and is therefore easier to use if you want different spatial
#' observations (and covariates) for each time slice. Eventually [sdmTMB_sim()]
#' will be depreciated in favour of this version.
#'
#' @param formula A *one-sided* formula describing the fixed-effect structure.
#' @param data A data frame containing the predictors described in `formula` and the
#'   time column if `time` is specified.
#' @param mesh Output from [make_mesh()].
#' @param time The time column name.
#' @param family Family as in [sdmTMB()].
#' @param B A vector of beta values (fixed-effect coefficient values).
#' @param range Parameter that controls the decay of spatial correlation. If a vector
#'  of length 2, `share_range` will be set to `FALSE` and the spatial and spatiotemporal
#'  ranges will be unique.
#' @param rho Spatiotemporal correlation between years; should be between -1 and
#'   1.
#' @param sigma_O SD of spatial process (Omega).
#' @param sigma_E SD of spatiotemporal process (Epsilon).
#' @param sigma_Z SD of spatially varying coefficient field (Zeta).
#' @param phi Observation error scale parameter (e.g., SD in Gaussian).
#' @param tweedie_p Tweedie p (power) parameter; between 1 and 2.
#' @param df Student-t degrees of freedom.
#' @param seed A value with which to set the random seed.
#' @param simulate_re Include random effect simulation?
#' @param previous_fit An optional previous [sdmTMB()] fit to pull parameter values.
#'   Will be over-ruled by any non-NULL specified parameter arguments.
#' @param ... Any other arguments to pass to [sdmTMB()].
#'
#' @return A data frame where:
#' * The 1st column is the time variable (if present).
#' * The 2nd and 3rd columns are the spatial coordinates.
#' * `omega_s` represents the simulated spatial random effects.
#' * `epsilon_st` represents the simulated spatiotemporal random effects.
#' * `eta` is the true value in link space
#' * `mu` is the true value in inverse link space.
#' * `observed` represents the simulated process with observation error.
#' * The remaining columns are the fixed-effect model matrix.
#' @export
#'
#' @examples
#' if (inla_installed()) {
#'   set.seed(123)
#'   # a1 is a fake predictor:
#'   predictor_dat <- data.frame(
#'     X = runif(300), Y = runif(300),
#'     a1 = rnorm(300), year = rep(1:6, each = 50)
#'   )
#'   mesh <- make_mesh(predictor_dat, xy_cols = c("X", "Y"), cutoff = 0.1)
#'
#'   sim_dat <- sdmTMB_simulate(
#'     formula = ~ 1 + a1,
#'     data = predictor_dat,
#'     time = "year",
#'     mesh = mesh,
#'     family = gaussian(link = "identity"),
#'     range = 0.5,
#'     sigma_E = 0.1,
#'     phi = 0.1,
#'     sigma_O = 0.2,
#'     seed = 3542,
#'     B = c(0.2, -0.4) # B0 = intercept, B1 = a1 slope
#'   )
#'
#'   fit <- sdmTMB(observed ~ a1, data = sim_dat, mesh = mesh, time = "year")
#'   fit
#'
#'   # example of supplying random field values:
#'   p <- predict(fit, newdata = NULL)
#'   sim_dat2 <- sdmTMB_simulate(
#'     formula = ~ 1 + a1,
#'     data = predictor_dat,
#'     time = "year",
#'     mesh = mesh,
#'     family = gaussian(link = "identity"),
#'     omega_s = p$omega_s,
#'     epsilon_st = p$epsilon_st,
#'     phi = 0.1,
#'     seed = 342,
#'     B = c(0.1, -0.2) # B0 = intercept, B1 = a1 slope
#'   )
#'
#'   if (require("ggplot2", quietly = TRUE)) {
#'     ggplot(sim_dat, aes(X, Y, colour = observed)) +
#'       geom_point() +
#'       facet_wrap(~year) +
#'       scale_color_gradient2()
#'   }
#' }
sdmTMB_simulate <- function(
  formula,
  data,
  mesh,
  family = gaussian(link = "identity"),
  time = NULL,
  previous_fit = NULL,
  B = NULL,
  range = NULL,
  rho = NULL,
  sigma_O = NULL,
  sigma_E = NULL,
  sigma_Z = NULL,
  phi = NULL,
  tweedie_p = NULL,
  df = NULL,
  simulate_re = TRUE,
  seed = sample.int(1e6, 1),
  ...) {

  if (!requireNamespace("INLA", quietly = TRUE)) {
    stop("INLA must be installed to use this function.", call. = FALSE)
  }

  if (is.null(previous_fit)) {
    assert_that(!is.null(range), !is.null(sigma_O) || !is.null(sigma_E), !is.null(B))
    if (!family$family %in% c("binomial", "poisson"))
      assert_that(!is.null(phi))
    assert_that(class(mesh) %in% "sdmTMBmesh")
    assert_that(tweedie_p > 1 && tweedie_p < 2 || is.null(tweedie_p))
    assert_that(df >= 1 || is.null(df))
    assert_that(all(range > 0))
    assert_that(length(range) %in% c(1, 2))
    assert_that(rho >= -1 && rho <= 1 || is.null(rho))
    assert_that(phi > 0 || is.null(phi))
    assert_that(sigma_O >= 0 || is.null(sigma_O))
    assert_that(all(sigma_E >= 0) || is.null(sigma_E))
    assert_that(sigma_Z >= 0 || is.null(sigma_Z))

    response <- get_response(formula)
    if (length(response) == 0L) {
      formula <- as.formula(paste("sdmTMB_response_", paste(as.character(formula), collapse = "")))
      data[["sdmTMB_response_"]] <- 0.1 # fake! does nothing but lets sdmTMB parse the formula
    }
  }

  # get tmb_data structure; parsed model matrices etc.:
  if (is.null(previous_fit)) {
    fit <- sdmTMB(
      formula = formula, data = data, mesh = mesh, time = time,
      family = family, do_fit = FALSE,
      share_range = length(range) == 1L,
      experimental = list(sim_re = simulate_re), ...
    )
  } else {
    fit <- previous_fit
  }
  params <- fit$tmb_params
  tmb_data <- fit$tmb_data

  if (!is.null(B)) {
    n_covariates <- length(B)
    assert_that(ncol(fit$tmb_data$X_ij) == length(B),
      msg = paste0("Number of specified fixed-effect `B` parameters does ",
        "not match model matrix columns implied by the formula.")
    )
  }

  if (is.null(sigma_O)) sigma_O <- 0
  if (is.null(sigma_Z)) sigma_Z <- 0
  if (is.null(sigma_E)) sigma_E <- 0

  if (!is.null(range)) {
    if (length(range) == 1L) range <- rep(range, 2)
    kappa <- sqrt(8)/range
    tau_O <- 1/(sqrt(4 * pi) * kappa[1] * sigma_O)
    tau_Z <- 1/(sqrt(4 * pi) * kappa[1] * sigma_Z)
    tau_E <- 1/(sqrt(4 * pi) * kappa[2] * sigma_E)
    params$ln_kappa <- log(kappa)
    params$ln_tau_O <- log(tau_O)
    params$ln_tau_E <- log(tau_E)
    params$ln_tau_Z <- log(tau_Z)
  }

  if (!is.null(B)) params$b_j <- B
  if (!is.null(phi)) params$ln_phi <- log(phi)
  if (!is.null(rho))  params$ar1_phi <- stats::qlogis(rho + 1)/2
  if (!is.null(df)) tmb_data$df <- df

  newobj <- TMB::MakeADFun(data = tmb_data, map = fit$tmb_map,
    random = fit$tmb_random, parameters = params, DLL = "sdmTMB")

  set.seed(seed)
  s <- newobj$simulate()

  d <- list()
  if (!is.null(time)) d[[time]] <- data[[time]]
  d[[mesh$xy_cols[1]]] <- data[[mesh$xy_cols[1]]]
  d[[mesh$xy_cols[2]]] <- data[[mesh$xy_cols[2]]]
  d[["omega_s"]] <- s$omega_s_A
  d[["epsilon_st"]] <- s$epsilon_st_A_vec
  d[["zeta_s"]] <- s$zeta_s_A
  d[["mu"]] <- s$eta_i
  d[["eta"]] <- family$linkfun(s$eta_i)
  d[["observed"]] <- s$y_i
  d <- do.call("data.frame", d)
  cbind(d, fit$tmb_data$X_ij)
}

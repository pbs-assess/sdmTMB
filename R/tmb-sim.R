#' Simulate from a spatial/spatiotemporal model
#'
#' @description
#' `r lifecycle::badge("experimental")`
#'
#' `sdmTMB_simulate()` (vs. [sdmTMB_sim()]) (1) uses TMB for simulation and is
#' therefore **much** faster and more flexible, (2) is set up to take a formula and
#' a data frame and is therefore easier to use if you want different spatial
#' observations (and covariates) for each time slice. [sdmTMB_sim()] uses pure
#' R for simulation. `sdmTMB_simulate()` uses TMB for simulation.
#'
#' [simulate.sdmTMB()], on the other hand, takes an existing fit and simulates
#' new observations and optionally new random effects.
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
#' @param fixed_re A list of optional random effects to fix at specified
#'    (e.g., previously estimated) values. Values of `NULL` will result
#'    in the random effects being simulated.
#' @param previous_fit (Depreciated; please use [simulate.sdmTMB()]).
#'   An optional previous [sdmTMB()] fit to pull parameter values.
#'   Will be over-ruled by any non-NULL specified parameter arguments.
#' @param seed Seed number.
#' @param ... Any other arguments to pass to [sdmTMB()].
#'
#' @return A data frame where:
#' * The 1st column is the time variable (if present).
#' * The 2nd and 3rd columns are the spatial coordinates.
#' * `omega_s` represents the simulated spatial random effects (only if present).
#' * `zeta_s` represents the simulated spatial varying covariate field (only if present).
#' * `epsilon_st` represents the simulated spatiotemporal random effects (only if present).
#' * `eta` is the true value in link space
#' * `mu` is the true value in inverse link space.
#' * `observed` represents the simulated process with observation error.
#' * The remaining columns are the fixed-effect model matrix.
#' @export
#' @seealso [simulate.sdmTMB()]
#'
#' @examples
#' if (inla_installed()) {
#'   set.seed(123)
#'
#'   # make fake predictor(s) (a1) and sampling locations:
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
#'     family = gaussian(),
#'     range = 0.5,
#'     sigma_E = 0.1,
#'     phi = 0.1,
#'     sigma_O = 0.2,
#'     seed = 42,
#'     B = c(0.2, -0.4) # B0 = intercept, B1 = a1 slope
#'   )
#'   head(sim_dat)
#'
#'   if (require("ggplot2", quietly = TRUE)) {
#'     ggplot(sim_dat, aes(X, Y, colour = observed)) +
#'       geom_point() +
#'       facet_wrap(~year) +
#'       scale_color_gradient2()
#'   }
#'
#'   # fit to the simulated data:
#'   fit <- sdmTMB(observed ~ a1, data = sim_dat, mesh = mesh, time = "year")
#'   fit
#'
#'   # example supplying previous fit, simulating new random effects,
#'   # and changing spatial SD (sigma_O) and observation error (phi):
#'   sim_dat2 <- sdmTMB_simulate(
#'     previous_fit = fit,
#'     simulate_re = TRUE, phi = 0.04, sigma_O = 0.4
#'   )
#'   head(sim_dat2)
#' }
sdmTMB_simulate <- function(formula,
                            data,
                            mesh,
                            family = gaussian(link = "identity"),
                            previous_fit = NULL,
                            time = NULL,
                            B = NULL,
                            range = NULL,
                            rho = NULL,
                            sigma_O = NULL,
                            sigma_E = NULL,
                            sigma_Z = NULL,
                            phi = NULL,
                            tweedie_p = NULL,
                            df = NULL,
                            fixed_re = list(omega_s = NULL, epsilon_st = NULL, zeta_s = NULL),
                            seed = sample.int(1e6, 1),
                            ...) {
  if (!requireNamespace("INLA", quietly = TRUE)) {
    nice_stop("INLA must be installed to use this function.")
  }

  if (!is.null(previous_fit)) mesh <- previous_fit$spde
  if (!is.null(previous_fit)) data <- previous_fit$data

  assert_that(tweedie_p > 1 && tweedie_p < 2 || is.null(tweedie_p))
  assert_that(df >= 1 || is.null(df))
  assert_that(all(range > 0) || is.null(range))
  assert_that(length(range) %in% c(1, 2) || is.null(range))
  assert_that(rho >= -1 && rho <= 1 || is.null(rho))
  assert_that(phi > 0 || is.null(phi))
  assert_that(sigma_O >= 0 || is.null(sigma_O))
  assert_that(all(sigma_E >= 0) || is.null(sigma_E))
  assert_that(all(sigma_Z >= 0) || is.null(sigma_Z))

  if (is.null(previous_fit)) {
    assert_that(is(mesh, "sdmTMBmesh"))
    assert_that(!is.null(range), !is.null(sigma_O) || !is.null(sigma_E), !is.null(B))
    if (!family$family %in% c("binomial", "poisson")) {
      assert_that(!is.null(phi))
    }

    response <- get_response(formula)
    if (length(response) == 0L) {
      formula <- as.formula(paste("sdmTMB_response_", paste(as.character(formula), collapse = "")))
      data[["sdmTMB_response_"]] <- 0.1 # fake! does nothing but lets sdmTMB parse the formula
      if (family$family %in% c("binomial", "poisson", "nbinom2", "nbinom1", "truncated_nbinom2", "truncated_nbinom1")) {
        data[["sdmTMB_response_"]] <- 1
      }
    }

    .sim_re <- list(
      omega = TRUE, epsilon = TRUE, zeta = TRUE,
      IID = TRUE, RW = TRUE, smooth = TRUE
    )
    if (!is.null(fixed_re$omega_s)) {
      .sim_re$omega <- FALSE
    }
    if (!is.null(fixed_re$epsilon_st)) {
      .sim_re$epsilon <- FALSE
    }
    if (!is.null(fixed_re$zeta_s)) {
      .sim_re$zeta <- FALSE
    }
    .sim_re <- as.integer(unlist(.sim_re))

    # get tmb_data structure; parsed model matrices etc.:
    fit <- sdmTMB(
      formula = formula, data = data, mesh = mesh, time = time,
      family = family, do_fit = FALSE,
      share_range = length(range) == 1L,
      experimental = list(sim_re = .sim_re), ...
    )
    params <- fit$tmb_params
  } else {
    fit <- previous_fit
    params <- fit$tmb_obj$env$parList()
  }
  tmb_data <- fit$tmb_data
  # tmb_data$sim_re <- as.integer(simulate_re)

  if (!is.null(B)) {
    n_covariates <- length(B)
    assert_that(ncol(fit$tmb_data$X_ij[[1]]) == length(B),
      msg = paste0(
        "Number of specified fixed-effect `B` parameters does ",
        "not match model matrix columns implied by the formula."
      )
    )
  }

  if (!is.null(previous_fit)) {
    range <- fit$tmb_obj$report()$range
  }

  if (is.null(previous_fit)) {
    if (is.null(sigma_O)) sigma_O <- 0
    if (is.null(sigma_Z)) sigma_Z <- matrix(0, nrow = 0L, ncol = 0L) # DELTA FIXME
    if (is.null(sigma_E)) sigma_E <- 0
  }

  if (length(range) == 1L) range <- rep(range, 2)

  kappa <- sqrt(8) / range
  params$ln_kappa <- matrix(log(kappa), ncol = 1L) # TODO DELTA

  if (!is.null(sigma_O) || is.null(previous_fit)) {
    tau_O <- 1 / (sqrt(4 * pi) * kappa[1] * sigma_O)
    params$ln_tau_O <- log(tau_O)
  }
  if (!is.null(sigma_Z) || is.null(previous_fit)) {
    tau_Z <- 1 / (sqrt(4 * pi) * kappa[1] * sigma_Z)
    params$ln_tau_Z <- matrix(log(tau_Z), nrow = nrow(sigma_Z), ncol = 1L) # DELTA FIXME
  }
  if (!is.null(sigma_E) || is.null(previous_fit)) {
    tau_E <- 1 / (sqrt(4 * pi) * kappa[2] * sigma_E)
    params$ln_tau_E <- log(tau_E)
  }

  if (!is.null(B)) params$b_j <- matrix(B, ncol = 1L) # TODO DELTA
  if (!is.null(phi)) params$ln_phi <- log(phi)
  if (!is.null(rho)) {
    if (rho != 0) tmb_data$ar1_fields <- 1L
    params$ar1_phi <- stats::qlogis((rho + 1) / 2)
  }
  if (!is.null(df)) tmb_data$df <- df
  if (!is.null(tweedie_p)) params$thetaf <- stats::qlogis(tweedie_p - 1)

  if (!is.null(fixed_re$omega_s)) {
    params$omega_s <- fixed_re$omega_s
  }
  if (!is.null(fixed_re$epsilon_st)) {
    params$epsilon_st <- fixed_re$epsilon_st
  }
  if (!is.null(fixed_re$zeta_s)) {
    params$zeta_s <- fixed_re$zeta_s
  }

  newobj <- TMB::MakeADFun(
    data = tmb_data, map = fit$tmb_map,
    random = fit$tmb_random, parameters = params, DLL = "sdmTMB"
  )

  set.seed(seed)
  s <- newobj$simulate()

  d <- list()
  if (!is.null(fit$time)) d[[fit$time]] <- data[[fit$time]]
  d[[mesh$xy_cols[1]]] <- data[[mesh$xy_cols[1]]]
  d[[mesh$xy_cols[2]]] <- data[[mesh$xy_cols[2]]]
  d[["omega_s"]] <- if (all(s$omega_s_A != 0)) s$omega_s_A
  d[["epsilon_st"]] <- if (all(s$epsilon_st_A_vec != 0)) s$epsilon_st_A_vec
  d[["zeta_s"]] <- if (all(s$zeta_s_A != 0)) s$zeta_s_A
  d[["mu"]] <- family$linkinv(s$eta_i)
  d[["eta"]] <- s$eta_i
  d[["observed"]] <- s$y_i
  d <- do.call("data.frame", d)
  cbind(d, fit$tmb_data$X_ij)
}

#' Simulate from a fitted sdmTMB model
#'
#' `simulate.sdmTMB` is an S3 method for producing a matrix of simulations from
#' a fitted model. This is similar to [lme4::simulate.merMod()] and
#' [glmmTMB::simulate.glmmTMB()]. It can be used with the \pkg{DHARMa} package
#' among other uses.
#'
#' @method simulate sdmTMB
#' @param object sdmTMB model
#' @param nsim Number of response lists to simulate. Defaults to 1.
#' @param seed Random number seed
#' @param params Whether the parameters used in the simulation should come from
#'   the Maximum Likelihood Estimate (`"mle"`) or from new draws from the joint
#'   precision matrix assuming they are multivariate normal distributed
#'   (`"mvn"`).
#' @param re_form `NULL` to specify a simulation conditional on fitted random
#'   effects (this only simulates observation error). `~0` or `NA` to simulate
#'   new random affects (smoothers, which internally are random effects, will
#'   not be simulated as new).
#' @param model If a delta/hurdle model, which model to simulate from?
#'   `NA` = combined, `1` = first model, `2` = second mdoel.
#' @param tmbstan_model An optional model fit via [tmbstan::tmbstan()]. If
#'   provided the parameters will be drawn from the MCMC samples and new
#'   observation error will be added. See the example in [extract_mcmc()].
#' @param ... Extra arguments (not used)
#' @return Returns a matrix; number of columns is `nsim`.
#' @importFrom stats simulate
#'
#' @seealso [sdmTMB_simulate()], [dharma_residuals()]
#'
#' @export
#' @examples
#' if (inla_installed()) {
#'
#' # start with some data simulated from scratch:
#' set.seed(1)
#' predictor_dat <- data.frame(X = runif(300), Y = runif(300), a1 = rnorm(300))
#' mesh <- make_mesh(predictor_dat, xy_cols = c("X", "Y"), cutoff = 0.1)
#' dat <- sdmTMB_simulate(
#'   formula = ~ 1 + a1,
#'   data = predictor_dat,
#'   mesh = mesh,
#'   family = poisson(),
#'   range = 0.5,
#'   sigma_O = 0.2,
#'   seed = 42,
#'   B = c(0.2, -0.4) # B0 = intercept, B1 = a1 slope
#' )
#' fit <- sdmTMB(observed ~ 1 + a1, data = dat, family = poisson(), mesh = mesh)
#'
#' # simulate from the model:
#' s1 <- simulate(fit, nsim = 300)
#' dim(s1)
#'
#' # test whether fitted models are consistent with the observed number of zeros:
#' sum(s1 == 0)/length(s1)
#' sum(dat$observed == 0) / length(dat$observed)
#'
#' # use the residuals with DHARMa:
#' if (require("DHARMa", quietly = TRUE)) {
#'   pred_fixed <- fit$family$linkinv(predict(fit)$est_non_rf)
#'   r <- DHARMa::createDHARMa(
#'     simulatedResponse = s1,
#'     observedResponse = dat$observed,
#'     fittedPredictedResponse = pred_fixed
#'   )
#'   plot(r)
#'   DHARMa::testResiduals(r)
#'   DHARMa::testSpatialAutocorrelation(r, x = dat$X, y = dat$Y)
#'   DHARMa::testZeroInflation(r)
#' }
#'
#' # simulate with the parameters drawn from the joint precision matrix:
#' s2 <- simulate(fit, nsim = 1, params = "MVN")
#'
#' # simulate with new random fields:
#' s3 <- simulate(fit, nsim = 1, re_form = ~ 0)
#'
#' # simulate with new random fields and new parameter draws:
#' s3 <- simulate(fit, nsim = 1, params = "MVN", re_form = ~ 0)
#'
#' # simulate from a Stan model fit with new observation error:
#' \donttest{
#' if (require("tmbstan", quietly = TRUE)) {
#'   stan_fit <- tmbstan::tmbstan(fit$tmb_obj, iter = 110, warmup = 100, chains = 1)
#'   # make sure `nsim` is <= number of samples from rstan
#'   s3 <- simulate(fit, nsim = 10, tmbstan_model = stan_fit)
#' }
#' }
#' }

simulate.sdmTMB <- function(object, nsim = 1L, seed = sample.int(1e6, 1L),
                            params = c("mle", "mvn"),
                            model = c(NA, 1, 2),
                            re_form = NULL, tmbstan_model = NULL, ...) {
  set.seed(seed)
  params <- tolower(params)
  params <- match.arg(params, choices = c("mle", "mvn"))

  # re_form stuff
  conditional_re <- !(!is.null(re_form) && ((re_form == ~0) || identical(re_form, NA)))
  tmb_dat <- object$tmb_data
  if (conditional_re) {
    tmb_dat$sim_re <- rep(0L, length(object$tmb_data$sim_re)) # don't simulate any REs
  } else {
    stopifnot(length(object$tmb_data$sim_re) == 6L) # in case this gets changed
    tmb_dat$sim_re <- c(rep(1L, 5L), 0L) # last is smoothers; don't simulate them
  }
  newobj <- TMB::MakeADFun(
    data = tmb_dat, map = object$tmb_map,
    random = object$tmb_random, parameters = object$tmb_obj$env$parList(), DLL = "sdmTMB"
  )

  # params MLE/MVN stuff
  if (is.null(tmbstan_model)) {
    if (params == "MVN") {
      new_par <- rmvnorm_prec(object$tmb_obj$env$last.par.best, object$sd_report, nsim)
    } else {
      new_par <- object$tmb_obj$env$last.par.best
    }
  } else {
    new_par <- extract_mcmc(tmbstan_model)
  }

  # do the sim
  if (params == "MVN" || !is.null(tmbstan_model)) { # we have a matrix
    ret <- lapply(seq_len(nsim), function(i) {
      newobj$simulate(par = new_par[, i, drop = TRUE], complete = FALSE)$y_i
    })
  } else {
    ret <- lapply(seq_len(nsim), function(i) newobj$simulate(par = new_par, complete = FALSE)$y_i)
  }

  if (isTRUE(object$family$delta)) {
    if (is.na(model[[1]])) {
      ret <- lapply(ret, function(.x) ifelse(!is.na(.x[,2]), .x[,2], .x[,1]))
    } else if (model[[1]] == 1) {
      ret <- lapply(ret, function(.x) .x[,1])
    } else if (model[[1]] == 2) {
      ret <- lapply(ret, function(.x) ifelse(!is.na(.x[,2]), .x[,2], NA))
    } else {
      abort("`model` argument isn't valid; should be NA, 1, or 2.")
    }
  }

  do.call(cbind, ret)
}

#' Get DHARMa residuals
#'
#' @param simulated_response Output from [simulate.sdmTMB()].
#' @param object Output from [sdmTMB()].
#' @param plot Logical.
#' @param ... Other arguments to pass to [DHARMa::createDHARMa()].
# @param fitted_column The column from the output of [predict.sdmTMB()] to pass
#   to [DHARMa::createDHARMa()]'s `fittedPredictedResponse` argument.
#'
#' @return
#' A data frame of observed and expected values is invisibly returned,
#' so you can set `plot = FALSE` and assign the output to an object if you wish
#' to plot the residuals yourself. See the examples.
#' @export
#'
#' @seealso [simulate.sdmTMB()], [residuals.sdmTMB()]
#'
#' @examples
#' if (inla_installed()) {
#' fit <- sdmTMB(density ~ as.factor(year) + s(depth, k = 3),
#'   data = pcod_2011, time = "year", mesh = pcod_mesh_2011,
#'   family = tweedie(link = "log"), spatial = "off",
#'   spatiotemporal = "off")
#'
#' # The `simulated_response` argument is first so the output from
#' # simulate() can be piped to dharma_residuals():
#' # simulate(fit, nsim = 500) %>% dharma_residuals(fit)
#'
#' s <- simulate(fit, nsim = 500)
#' dharma_residuals(s, fit)
#' r <- dharma_residuals(s, fit, plot = FALSE)
#' head(r)
#' plot(r$expected, r$observed)
#' abline(a = 0, b = 1)
#' }

dharma_residuals <- function(simulated_response, object, plot = TRUE, ...) {
  if (!requireNamespace("DHARMa", quietly = TRUE)) {
    nice_stop("DHARMa must be installed to use this function.")
  }

  assert_that(class(object) == "sdmTMB")
  assert_that(is.logical(plot))
  assert_that(is.matrix(simulated_response))
  assert_that(nrow(simulated_response) == nrow(object$response))
  if (isTRUE(object$family$delta)) {
    y <- ifelse(!is.na(object$response[,2]),
      object$response[,2], object$response[,1])
  } else {
    y <- object$response[,1]
  }
  y <- as.numeric(y)

  n_orig <- suppressWarnings(TMB::openmp(NULL))
  if (n_orig > 0 && .Platform$OS.type == "unix") { # openMP is supported
    TMB::openmp(n = 1L)
    on.exit({TMB::openmp(n = n_orig)})
  }

  p <- predict(object, type = "response")
  # fitted <- object$family$linkinv(p[["est_non_rf"]])
  fitted <- p$est
  res <- DHARMa::createDHARMa(
    simulatedResponse = simulated_response,
    observedResponse = y,
    fittedPredictedResponse = fitted,
    ...
  )
  u <- res$scaledResiduals
  n <- length(u)
  m <- seq_len(n) / (n + 1)
  z <- stats::qqplot(m, u, plot.it = FALSE)
  if (plot) {
    DHARMa::plotQQunif(
      res,
      testUniformity = FALSE,
      testOutliers = FALSE, testDispersion = FALSE
    )
  }
  invisible(data.frame(observed = z$y, expected = z$x))
}


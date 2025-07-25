#' Simulate from a spatial/spatiotemporal model
#'
#' `sdmTMB_simulate()` uses TMB to simulate *new* data given specified parameter
#' values. [simulate.sdmTMB()], on the other hand, takes an *existing* model fit
#' and simulates new observations and optionally new random effects.
#'
#' @param formula A *one-sided* formula describing the fixed-effect structure.
#'   Random intercepts are not (yet) supported. Fixed effects should match
#'   the corresponding `B` argument vector of coefficient values.
#' @param data A data frame containing the predictors described in `formula` and
#'   the time column if `time` is specified.
#' @param mesh Output from [make_mesh()].
#' @param time The time column name.
#' @param family Family as in [sdmTMB()]. Delta families are not supported.
#'   Instead, simulate the two component models separately and combine.
#' @param B A vector of beta values (fixed-effect coefficient values).
#' @param range Parameter that controls the decay of spatial correlation. If a
#'   vector of length 2, `share_range` will be set to `FALSE` and the spatial
#'   and spatiotemporal ranges will be unique.
#' @param rho Spatiotemporal correlation between years; should be between -1 and
#'   1.
#' @param sigma_O SD of spatial process (Omega).
#' @param sigma_E SD of spatiotemporal process (Epsilon).
#' @param sigma_Z SD of spatially varying coefficient field (Zeta).
#' @param phi Observation error scale parameter (e.g., SD in Gaussian).
#' @param tweedie_p Tweedie p (power) parameter; between 1 and 2.
#' @param df Student-t degrees of freedom.
#' @param threshold_coefs An optional vector of threshold coefficient values
#'   if the `formula` includes `breakpt()` or `logistic()`. If `breakpt()`,
#'   these are slope and cut values. If `logistic()`, these are the threshold at
#'   which the function is 50% of the maximum, the threshold at which the
#'   function is 95% of the maximum, and the maximum. See the model description
#'   vignette for details.
#' @param fixed_re A list of optional random effects to fix at specified
#'    (e.g., previously estimated) values. Values of `NULL` will result
#'    in the random effects being simulated.
#' @param previous_fit (**Deprecated**; please use [simulate.sdmTMB()]).
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
#
#   # example supplying previous fit, simulating new random effects,
#   # and changing spatial SD (sigma_O) and observation error (phi):
#   sim_dat2 <- sdmTMB_simulate(
#     previous_fit = fit,
#     simulate_re = TRUE, phi = 0.04, sigma_O = 0.4
#   )
#   head(sim_dat2)
sdmTMB_simulate <- function(formula,
                            data,
                            mesh,
                            family = gaussian(link = "identity"),
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
                            threshold_coefs = NULL,
                            fixed_re = list(omega_s = NULL, epsilon_st = NULL, zeta_s = NULL),
                            previous_fit = NULL,
                            seed = sample.int(1e6, 1),
                            ...) {

  if (!is.null(previous_fit)) stop("`previous_fit` is deprecated. See `simulate.sdmTMB()`", call. = FALSE)
  if (!is.null(previous_fit)) mesh <- previous_fit$spde
  if (!is.null(previous_fit)) data <- previous_fit$data
  # if (!missing(seed)) {
  #   msg <- c("The `seed` argument may be deprecated in the future.",
  #     "We recommend instead setting the seed manually with `set.seed()` prior to calling `sdmTMB_simulate()`.",
  #     "We have encountered some situations where setting the seed via this argument does not have the intended effect.")
  #   cli_inform(msg)
  # }

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
      # experimental = list(sim_re = .sim_re),
      ...
    )
    params <- fit$tmb_params
  } else {
    fit <- previous_fit
    params <- fit$tmb_obj$env$parList()
  }
  tmb_data <- fit$tmb_data
  tmb_data$sim_re <- as.integer(.sim_re)
  # tmb_data$sim_re <- c(1L, 0L, 0L, 0L, 0L, 0L)

  if (!is.null(B)) {
    n_covariates <- length(B)
    assert_that(ncol(fit$tmb_data$X_ij[[1]]) == length(B),
      msg = paste0(
        "Number of specified fixed-effect `B` parameters does ",
        "not match model matrix columns implied by the formula."
      )
    )
  }

  if (tmb_data$threshold_func > 0) {
    if (is.null(threshold_coefs)) {
      cli::cli_abort("Break point or logistic formula detected without `threshold_coefs` defined.")
    }
  }
  if (!is.null(threshold_coefs)) {
    if(!is.matrix(threshold_coefs)) threshold_coefs <- matrix(threshold_coefs, ncol=1)
    params$b_threshold <- threshold_coefs
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
    if (rho != 0 && rho < 1) {
      tmb_data$ar1_fields <- 1L
      params$ar1_phi <- stats::qlogis((rho + 1) / 2)
    } else if (rho == 1) {
      tmb_data$rw_fields <- 1L
    }
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
    random = fit$tmb_random, parameters = params, DLL = "sdmTMB",
    checkParameterOrder = FALSE
  )

  set.seed(seed)
  s <- newobj$simulate()

  d <- list()
  if (!is.null(fit$time)) d[[fit$time]] <- data[[fit$time]]
  d[[mesh$xy_cols[1]]] <- data[[mesh$xy_cols[1]]]
  d[[mesh$xy_cols[2]]] <- data[[mesh$xy_cols[2]]]

  d[["omega_s"]] <- if (sum(sigma_O) > 0) s$omega_s_A
  d[["epsilon_st"]] <- if (sum(sigma_E) > 0) s$epsilon_st_A_vec
  d[["zeta_s"]] <- if (sum(sigma_Z) > 0) s$zeta_s_A

  # # Warnings for fields collapsing to 0
  # info_collapse <- function(sig, vec, .par, .name) {
  #   if (sum(sig) > 0 && all (vec == 0)) {
  #     msg <- paste0("The ", .name, " has been returned as all zeros although ", .par,
  #       " was specified as > 0. Try making your mesh finer, e.g., with a lower ",
  #       "`cutoff` or a higher number of knots. Triangle edge length needs to be ",
  #       "lower than the range size (distance correlation is effectively independent.")
  #     cli::cli_alert_info(msg)
  #   }
  # }
  # info_collapse(sigma_O, s$omega_s_A, "sigma_O", "spatial field")
  # info_collapse(sigma_E, s$epsilon_st_A_vec, "sigma_E", "spatiotemporal field")
  # info_collapse(sigma_Z, s$zeta_s_A, "sigma_Z", "spatially varying coefficient field")

  if (any(family$family %in% c("truncated_nbinom1", "truncated_nbinom2"))) {
    d[["mu"]] <- family$linkinv(s$eta_i, phi = phi)
  } else {
    d[["mu"]] <- family$linkinv(s$eta_i)
  }
  d[["eta"]] <- s$eta_i
  d[["observed"]] <- s$y_i
  d <- do.call("data.frame", d)
  d <- cbind(d, fit$tmb_data$X_ij)

  tpar <- fit$threshold_parameter
  if (tmb_data$threshold_func == 1L) {
    d[[paste0(tpar, "-slope")]] <- threshold_coefs[[1]]
    d[[paste0(tpar, "-breakpt")]] <- threshold_coefs[[2]]
  }
  if (tmb_data$threshold_func == 2L) {
    d[[paste0(tpar, "-s50")]] <- threshold_coefs[[1]]
    d[[paste0(tpar, "-s95")]] <- threshold_coefs[[2]]
    d[[paste0(tpar, "-smax")]] <- threshold_coefs[[3]]
  }

  d
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
#' @param type How parameters should be treated. `"mle-eb"`: fixed effects
#'   are at their maximum likelihood (MLE) estimates  and random effects are at
#'   their empirical Bayes (EB) estimates. `"mle-mvn"`: fixed effects are at
#'   their MLEs but random effects are taken from a single approximate sample.
#'   This latter option is a suggested approach if these simulations will be
#'   used for goodness of fit testing (e.g., with the DHARMa package).
#' @param re_form `NULL` to specify a simulation conditional on fitted random
#'   effects (this only simulates observation error). `~0` or `NA` to simulate
#'   new random affects (smoothers, which internally are random effects, will
#'   not be simulated as new).
#' @param mle_mvn_samples Applies if `type = "mle-mvn"`. If `"single"`, take
#'   a single MVN draw from the random effects. If `"multiple"`, take an MVN
#'   draw from the random effects for each of the `nsim`.
#' @param model If a delta/hurdle model, which model to simulate from?
#'   `NA` = combined, `1` = first model, `2` = second mdoel.
#' @param newdata Optional new data frame from which to simulate.
#' @param mcmc_samples An optional matrix of MCMC samples. See `extract_mcmc()`
#'   in the \href{https://github.com/pbs-assess/sdmTMBextra}{sdmTMBextra}
#'   package.
#' @param return_tmb_report Return the \pkg{TMB} report from `simulate()`? This
#'   lets you parse out whatever elements you want from the simulation.
#'   Not usually needed.
#' @param observation_error Logical. Simulate observation error?
#' @param size A vector of size (trials) in the case of a binomial family with
#'   `newdata` specified. If left `NULL` with `newdata`, will be assumed to
#'   be a vector of 1s.
#' @param silent Logical. Silent?
#' @param ... Extra arguments passed to [predict.sdmTMB()]. E.g., one may wish
#'   to pass an `offset` argument if `newdata` are supplied in a model with an
#'   offset.
#' @return Returns a matrix; number of columns is `nsim`.
#' @importFrom stats simulate
#'
#' @seealso [sdmTMB_simulate()]
#'
#' @export
#' @examples
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
#' # simulate with random effects sampled from their approximate posterior
#' s2 <- simulate(fit, nsim = 1, params = "mle-mvn")
#' # these may be useful in conjunction with DHARMa simulation-based residuals
#'
#' # simulate with new random fields:
#' s3 <- simulate(fit, nsim = 1, re_form = ~ 0)

simulate.sdmTMB <- function(object, nsim = 1L, seed = sample.int(1e6, 1L),
                            type = c("mle-eb", "mle-mvn"),
                            model = c(NA, 1, 2),
                            newdata = NULL,
                            re_form = NULL,
                            mle_mvn_samples = c("single", "multiple"),
                            mcmc_samples = NULL,
                            return_tmb_report = FALSE,
                            observation_error = TRUE,
                            size = NULL,
                            silent = FALSE,
                            ...) {
  set.seed(seed)
  type <- tolower(type)
  type <- match.arg(type)
  mle_mvn_samples <- match.arg(mle_mvn_samples)
  assert_that(as.integer(model[[1]]) %in% c(NA_integer_, 1L, 2L))

  # need to re-attach environment if in fresh session
  reinitialize(object)

  if (is.null(object$tmb_random) && type == "mle-mvn") {
    type <- "mle-eb" # no random effects to sample from
  }

  # re_form stuff
  conditional_re <- !(!is.null(re_form) && ((re_form == ~0) || identical(re_form, NA)))
  tmb_dat <- object$tmb_data
  if (conditional_re) {
    tmb_dat$sim_re <- rep(0L, length(object$tmb_data$sim_re)) # don't simulate any REs
  } else {
    stopifnot(length(object$tmb_data$sim_re) == 6L) # in case this gets changed
    tmb_dat$sim_re <- c(rep(1L, 5L), 0L) # last is smoothers; don't simulate them
  }

  if (!is.null(newdata)) {
    # generate prediction TMB data list
    p <- predict(object, newdata = newdata, return_tmb_data = TRUE, ...)
    # move data elements over
    p <- move_proj_to_tmbdat(p, object, newdata, called_by_simulate = TRUE, size = size)
    p$sim_re <- tmb_dat$sim_re
    tmb_dat <- p
  }

  tmb_dat$sim_obs <- as.integer(observation_error)

  newobj <- TMB::MakeADFun(
    data = tmb_dat, map = object$tmb_map,
    random = object$tmb_random, parameters = object$tmb_obj$env$parList(), DLL = "sdmTMB"
  )

  # params MLE/MVN stuff
  if (is.null(mcmc_samples)) {
    if (type == "mle-mvn") {
      if (mle_mvn_samples == "single") {
        new_par <- .one_sample_posterior(object)
        new_par <- replicate(nsim, new_par)
      } else {
        new_par <- lapply(seq_len(nsim), \(i) .one_sample_posterior(object))
        new_par <- do.call(cbind, new_par)
      }
    } else if (type == "mle-eb") {
      new_par <- object$tmb_obj$env$last.par.best
      new_par <- lapply(seq_len(nsim), \(i) new_par)
      new_par <- do.call(cbind, new_par)
    } else {
      cli_abort("`type` type not defined")
    }
  } else {
    new_par <- mcmc_samples
  }

  # do the simulation
  if (!silent) cli::cli_progress_bar("Simulating", total = nsim)
  ret <- list()
  if (!is.null(mcmc_samples)) { # we have a matrix
    for (i in seq_len(nsim)) {
      if (!silent) cli::cli_progress_update()
      ret[[i]] <- newobj$simulate(par = new_par[, i, drop = TRUE], complete = FALSE)
      if (!return_tmb_report) ret[[i]] <- ret[[i]]$y_i
    }
  } else {
    for (i in seq_len(nsim)) {
      if (!silent) cli::cli_progress_update()
      ret[[i]] <- newobj$simulate(par = new_par[, i, drop = TRUE], complete = FALSE)
      if (!return_tmb_report) ret[[i]] <- ret[[i]]$y_i
    }
  }
  if (!silent) cli::cli_progress_done()

  if (!return_tmb_report) {
    if (is_delta(object)) {
      if (object$family[[1]]$family == "binomial" &&
          object$family$type == "standard" &&
          !observation_error) {
        ret <- lapply(ret, function(.x) {
          .x[,1] <- stats::plogis(.x[,1]) # using robust dbinom in logit space
          .x
        })
      }
      if (is.na(model[[1]])) {
        ret <- lapply(ret, function(.x) .x[,1] * .x[,2])
      } else if (model[[1]] == 1) {
        ret <- lapply(ret, function(.x) .x[,1])
      } else if (model[[1]] == 2) {
        ret <- lapply(ret, function(.x) .x[,2])
      } else {
        cli_abort("`model` argument isn't valid; should be NA, 1, or 2.")
      }
    }

    if (!is_delta(object) && !observation_error) {
      if (object$family$family == "binomial") {
        # using robust dbinom in logit space
        ret <- lapply(ret, function(.x) stats::plogis(.x))
      }
    }

    ret <- do.call(cbind, ret)

    if (!is.null(newdata)) {
      rownames(ret) <- newdata[[object$time]] # for use in index calcs
      attr(ret, "time") <- object$time
      if (is_delta(object)) {
        attr(ret, "link") <- object$family[[2]]$link
      } else {
        attr(ret, "link") <- object$family$link
      }
    }

    attr(ret, "type") <- type
    attr(ret, "mle_mvn_samples") <- mle_mvn_samples
  }
  ret
}

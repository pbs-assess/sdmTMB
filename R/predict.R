# @importFrom stats predict
# @rdname predict
# @export

#' Predict from an sdmTMB model
#'
#' Make predictions from an sdmTMB model; can predict on the original or new
#' data.
#'
#' @param object An object from [sdmTMB()].
#' @param newdata A data frame to make predictions on. This should be a data
#'   frame with the same predictor columns as in the fitted data and a time
#'   column (if this is a spatiotemporal model) with the same name as in the
#'   fitted data. There should be predictor data for each year in the original
#'   data set.
#' @param se_fit Should standard errors on predictions at the new locations
#'   given by `newdata` be calculated? Warning: the current implementation can
#'   be very slow for large data sets or high-resolution projections. A *much*
#'   faster option is to use the `nsim` argument below and calculate uncertainty
#'   on the simulations from the joint precision matrix.
#' @param return_tmb_object Logical. If `TRUE`, will include the TMB object in a
#'   list format output. Necessary for the [get_index()] or [get_cog()]
#'   functions.
#' @param area A vector of areas for survey grid cells. Only necessary if the
#'   output will be passed to [get_index()] or [get_cog()] and not all grid
#'   cells are of area 1. Should be the same length as the number of rows of
#'   `newdata`. If length 1, will be repeated to match the rows of data.
#' @param re_form `NULL` to specify including all spatial/spatiotemporal random
#'   effects in predictions. `~0` or `NA` for population-level predictions. Note
#'   that unlike lme4 or glmmTMB, this only affects what the standard errors are
#'   calculated on if `se_fit = TRUE`. This does not affect [get_index()]
#'   calculations.
#' @param re_form_iid `NULL` to specify including all random intercepts in the
#'   predictions. `~0` or `NA` for population-level predictions. No other
#'   options (e.g., some but not all random intercepts) are implemented yet.
#'   Only affects predictions with `newdata`. This also affects [get_index()].
#' @param nsim **Experimental.** If > 0, simulate from the joint precision matrix with `sims`
#'   draws Returns a matrix of `nrow(data)` by `sim` representing the estimates
#'   of the linear predictor (i.e., in link space). Can be useful for deriving
#'   uncertainty on predictions (e.g., `apply(x, 1, sd)`) or propagating
#'   uncertainty. This is currently the fastest way to generate estimates of
#'   uncertainty on predictions in space with sdmTMB.
#' @param sims **Deprecated**. Please use `nsim` instead.
#' @param sims_var **Experimental.** Which TMB reported variable from the model
#'   should be extracted from the joint precision matrix simulation draws?
#'   Defaults to the link-space predictions. Options include: `"omega_s"`,
#'   `"zeta_s"`, `"epsilon_st"`, and `"est_rf"` (as described below).
#'   Other options will be passed verbatim.
#' @param tmbstan_model A model fit with [tmbstan::tmbstan()]. See
#'   [extract_mcmc()] for more details and an example. If specified, the
#'   predict function will return a matrix of a similar form as if `nsim > 0`
#'   but representing Bayesian posterior samples from the Stan model.
#' @param ... Not implemented.
#'
#' @return
#' If `return_tmb_object = FALSE` (and `nsim = 0` and `tmbstan_model = NULL`):
#'
#' A data frame:
#' * `est`: Estimate in link space (everything is in link space)
#' * `est_non_rf`: Estimate from everything that isn't a random field
#' * `est_rf`: Estimate from all random fields combined
#' * `omega_s`: Spatial (intercept) random field that is constant through time
#' * `zeta_s`: Spatial slope random field
#' * `epsilon_st`: Spatiotemporal (intercept) random fields, could be
#'    off (zero), IID, AR1, or random walk
#'
#' If `return_tmb_object = TRUE` (and `nsim = 0` and `tmbstan_model = NULL`):
#'
#' A list:
#' * `data`: The data frame described above
#' * `report`: The TMB report on parameter values
#' * `obj`: The TMB object returned from the prediction run
#' * `fit_obj`: The original TMB model object
#'
#' In this case, you likely only need the `data` element as an end user.
#' The other elements are included for other functions.
#'
#' If `nsim > 0` or `tmbstan_model` is not `NULL`:
#'
#' A matrix:
#'
#' * Columns represent samples
#' * Rows represent predictions with one row per row of `newdata`
#'
#' @export
#'
#' @examples
#' if (require("ggplot2", quietly = TRUE) && inla_installed()) {
#'
#' d <- pcod_2011
#' mesh <- make_mesh(d, c("X", "Y"), cutoff = 30) # a coarse mesh for example speed
#' m <- sdmTMB(
#'  data = d, formula = density ~ 0 + as.factor(year) + depth_scaled + depth_scaled2,
#'  time = "year", mesh = mesh, family = tweedie(link = "log")
#' )
#'
#' # Predictions at original data locations -------------------------------
#'
#' predictions <- predict(m)
#' head(predictions)
#'
#' predictions$resids <- residuals(m, type = "randomized-quantile")
#'
#' ggplot(predictions, aes(X, Y, col = resids)) + scale_colour_gradient2() +
#'   geom_point() + facet_wrap(~year)
#' hist(predictions$resids)
#' qqnorm(predictions$resids);abline(a = 0, b = 1)
#'
#' # Predictions onto new data --------------------------------------------
#'
#' qcs_grid_2011 <- subset(qcs_grid, year >= min(pcod_2011$year))
#' predictions <- predict(m, newdata = qcs_grid_2011)
#'
#' # A short function for plotting our predictions:
#' plot_map <- function(dat, column = "est") {
#'   ggplot(dat, aes_string("X", "Y", fill = column)) +
#'     geom_raster() +
#'     facet_wrap(~year) +
#'     coord_fixed()
#' }
#'
#' plot_map(predictions, "exp(est)") +
#'   scale_fill_viridis_c(trans = "sqrt") +
#'   ggtitle("Prediction (fixed effects + all random effects)")
#'
#' plot_map(predictions, "exp(est_non_rf)") +
#'   ggtitle("Prediction (fixed effects and any time-varying effects)") +
#'   scale_fill_viridis_c(trans = "sqrt")
#'
#' plot_map(predictions, "est_rf") +
#'   ggtitle("All random field estimates") +
#'   scale_fill_gradient2()
#'
#' plot_map(predictions, "omega_s") +
#'   ggtitle("Spatial random effects only") +
#'   scale_fill_gradient2()
#'
#' plot_map(predictions, "epsilon_st") +
#'   ggtitle("Spatiotemporal random effects only") +
#'   scale_fill_gradient2()
#'
#' # Visualizing a marginal effect ----------------------------------------
#'
#' nd <- data.frame(depth_scaled =
#'   seq(min(d$depth_scaled), max(d$depth_scaled), length.out = 100))
#' nd$depth_scaled2 <- nd$depth_scaled^2
#'
#' # Because this is a spatiotemporal model, you'll need at least one time
#' # element. If time isn't also a fixed effect then it doesn't matter what you pick:
#' nd$year <- 2011L # L: integer to match original data
#' p <- predict(m, newdata = nd, se_fit = TRUE, re_form = NA)
#' ggplot(p, aes(depth_scaled, exp(est),
#'   ymin = exp(est - 1.96 * est_se), ymax = exp(est + 1.96 * est_se))) +
#'   geom_line() + geom_ribbon(alpha = 0.4)
#'
#' # Plotting marginal effect of a spline ---------------------------------
#'
#' m_gam <- sdmTMB(
#'  data = d, formula = density ~ 0 + as.factor(year) + s(depth_scaled, k = 5),
#'  time = "year", mesh = mesh, family = tweedie(link = "log")
#' )
#' plot_smooth(m_gam)
#'
#' # or manually:
#' nd <- data.frame(depth_scaled =
#'   seq(min(d$depth_scaled), max(d$depth_scaled), length.out = 100))
#' nd$year <- 2011L
#' p <- predict(m_gam, newdata = nd, se_fit = TRUE, re_form = NA)
#' ggplot(p, aes(depth_scaled, exp(est),
#'   ymin = exp(est - 1.96 * est_se), ymax = exp(est + 1.96 * est_se))) +
#'   geom_line() + geom_ribbon(alpha = 0.4)
#'
#' # Forecasting ----------------------------------------------------------
#' mesh <- make_mesh(d, c("X", "Y"), cutoff = 15)
#'
#' unique(d$year)
#' m <- sdmTMB(
#'   data = d, formula = density ~ 1,
#'   spatiotemporal = "AR1", # using an AR1 to have something to forecast with
#'   extra_time = 2019L, # `L` for integer to match our data
#'   spatial = "off",
#'   time = "year", mesh = mesh, family = tweedie(link = "log")
#' )
#'
#' # Add a year to our grid:
#' grid2019 <- qcs_grid_2011[qcs_grid_2011$year == max(qcs_grid_2011$year), ]
#' grid2019$year <- 2019L # `L` because `year` is an integer in the data
#' qcsgrid_forecast <- rbind(qcs_grid_2011, grid2019)
#'
#' predictions <- predict(m, newdata = qcsgrid_forecast)
#' plot_map(predictions, "exp(est)") +
#'   scale_fill_viridis_c(trans = "log10")
#' plot_map(predictions, "epsilon_st") +
#'   scale_fill_gradient2()
#'
#' # Estimating local trends ----------------------------------------------
#'
#' d <- pcod
#' d$year_scaled <- as.numeric(scale(d$year))
#' mesh <- make_mesh(pcod, c("X", "Y"), cutoff = 25)
#' m <- sdmTMB(data = d, formula = density ~ depth_scaled + depth_scaled2,
#'   mesh = mesh, family = tweedie(link = "log"),
#'   spatial_varying = ~ 0 + year_scaled, time = "year", spatiotemporal = "off")
#' nd <- qcs_grid
#' nd$year_scaled <- (nd$year - mean(d$year)) / sd(d$year)
#' p <- predict(m, newdata = nd)
#'
#' plot_map(p, "zeta_s") +
#'   ggtitle("Spatial slopes") +
#'   scale_fill_gradient2()
#'
#' plot_map(p, "est_rf") +
#'   ggtitle("Random field estimates") +
#'   scale_fill_gradient2()
#'
#' plot_map(p, "exp(est_non_rf)") +
#'   ggtitle("Prediction (fixed effects only)") +
#'   scale_fill_viridis_c(trans = "sqrt")
#'
#' plot_map(p, "exp(est)") +
#'   ggtitle("Prediction (fixed effects + all random effects)") +
#'   scale_fill_viridis_c(trans = "sqrt")
#' }

predict.sdmTMB <- function(object, newdata = object$data, se_fit = FALSE,
  return_tmb_object = FALSE,
  area = 1, re_form = NULL, re_form_iid = NULL, nsim = 0,
  sims = deprecated(), sims_var = "est", tmbstan_model = NULL, ...) {

  # check if this is visreg:
  sys_calls <- unlist(lapply(sys.calls(), deparse)) # retrieve function that called this
  visreg_df <- FALSE
  if (!is.null(visreg_df) && any(substr(sys_calls, 1, 6) == "visreg")) {
    visreg_df <- TRUE
    if (any(sys_calls == "residuals(fit)")) visreg_df <- FALSE
    # turn on standard error if in a function call
    indx <- which(substr(sys_calls,1,10) == "visregPred")
    if (length(indx) > 0 && any(unlist(strsplit(sys_calls[indx], ",")) == " se.fit = TRUE")) se_fit <- TRUE
  }

  if ("version" %in% names(object)) {
    check_sdmTMB_version(object$version)
  } else {
    stop("This looks like a very old version of a model fit. Re-fit the model before predicting with it.",
      call. = FALSE)
  }
  if (!"xy_cols" %in% names(object$spde)) {
    warning("It looks like this model was fit with make_spde(). ",
    "Using `xy_cols`, but future versions of sdmTMB may not be compatible with this.",
    "Please replace make_spde() with make_mesh().", call. = FALSE)
  } else {
    xy_cols <- object$spde$xy_cols
  }

  if (is_present(sims)) {
    deprecate_warn("0.0.21", "predict.sdmTMB(sims)", "predict.sdmTMB(nsim)")
  } else {
    sims <- nsim
  }

  n_orig <- suppressWarnings(TMB::openmp(NULL))
  if (n_orig > 0 && .Platform$OS.type == "unix") { # openMP is supported
    TMB::openmp(n = object$control$parallel)
    on.exit({TMB::openmp(n = n_orig)})
  }

  # from glmmTMB:
  pop_pred <- (!is.null(re_form) && ((re_form == ~0) || identical(re_form, NA)))
  pop_pred_iid <- (!is.null(re_form_iid) && ((re_form_iid == ~0) || identical(re_form_iid, NA)))
  if (pop_pred_iid) {
    exclude_RE <- rep(1L, length(object$tmb_data$exclude_RE))
  } else {
    exclude_RE <- object$tmb_data$exclude_RE
  }

  tmb_data <- object$tmb_data
  tmb_data$do_predict <- 1L

  if (!is.null(newdata)) {
    if (any(!xy_cols %in% names(newdata)) && isFALSE(pop_pred))
      stop("`xy_cols` (the column names for the x and y coordinates) ",
        "are not in `newdata`. Did you miss specifying the argument ",
        "`xy_cols` to match your data? The newer `make_mesh()` ",
        "(vs. `make_spde()`) takes care of this for you.", call. = FALSE
      )

    if (object$time == "_sdmTMB_time") newdata[[object$time]] <- 0L

    check_time_class(object, newdata)
    original_time <- as.numeric(sort(unique(object$data[[object$time]])))
    new_data_time <- as.numeric(sort(unique(newdata[[object$time]])))

    if (!all(new_data_time %in% original_time))
      stop("Some new time elements were found in `newdata`. ",
        "For now, make sure only time elements from the original dataset ",
        "are present. If you would like to predict on new time elements, see ",
        "the `extra_time` argument in `?predict.sdmTMB`.",
        call. = FALSE
      )

    if (!identical(new_data_time, original_time) & isFALSE(pop_pred)) {
      stop("The time elements in `newdata` are not identical to those ",
        "in the original dataset. For now, please predict on all time ",
        "elements and filter out those you don't need after. Please ",
        "let us know on the GitHub issues if this is important to you.",
        call. = FALSE
      )
    }

    # If making population predictions (with standard errors), we don't need
    # to worry about space, so fill in dummy values if the user hasn't made any:
    fake_spatial_added <- FALSE
    if (pop_pred) {
      for (i in c(1, 2)) {
        if (!xy_cols[[i]] %in% names(newdata)) {
          newdata[[xy_cols[[i]]]] <- mean(object$data[[xy_cols[[i]]]], na.rm = TRUE)
          fake_spatial_added <- TRUE
        }
      }
    }

    if (sum(is.na(new_data_time)) > 0)
      stop("There is at least one NA value in the time column. ",
        "Please remove it.", call. = FALSE)

    newdata$sdm_orig_id <- seq(1, nrow(newdata))
    fake_newdata <- unique(newdata[,xy_cols])
    fake_newdata[["sdm_spatial_id"]] <- seq(1, nrow(fake_newdata)) - 1L

    newdata <- base::merge(newdata, fake_newdata, by = xy_cols,
      all.x = TRUE, all.y = FALSE)
    newdata <- newdata[order(newdata$sdm_orig_id),, drop = FALSE]


    proj_mesh <- INLA::inla.spde.make.A(object$spde$mesh,
      loc = as.matrix(fake_newdata[,xy_cols, drop = FALSE]))

    # this formula has breakpt() etc. in it:
    thresh <- check_and_parse_thresh_params(object$formula, newdata)
    formula <- thresh$formula # this one does not

    nd <- newdata
    response <- get_response(object$formula)
    sdmTMB_fake_response <- FALSE
    if (!response %in% names(nd)) {
      nd[[response]] <- 0 # fake for model.matrix
      sdmTMB_fake_response <- TRUE
    }

    # deal with prediction IID random intercepts:
    RE_names <- barnames(object$split_formula$reTrmFormulas)
    ## not checking so that not all factors need to be in prediction:
    # fct_check <- vapply(RE_names, function(x) check_valid_factor_levels(data[[x]], .name = x), TRUE)
    proj_RE_indexes <- vapply(RE_names, function(x) as.integer(nd[[x]]) - 1L, rep(1L, nrow(nd)))

    if (!"mgcv" %in% names(object)) object[["mgcv"]] <- FALSE
    f2 <- remove_s_and_t2(object$split_formula$fixedFormula)
    tt <- stats::terms(f2)
    Terms <- stats::delete.response(tt)
    mf <- model.frame(Terms, newdata, xlev = object$xlevels)
    proj_X_ij <- model.matrix(Terms, mf, contrasts.arg = object$contrasts)

    sm <- parse_smoothers(object$formula, data = object$data, newdata = nd)

    if (!is.null(object$time_varying))
      proj_X_rw_ik <- model.matrix(object$time_varying, data = nd)
    else
      proj_X_rw_ik <- matrix(0, ncol = 1, nrow = 1) # dummy

    if (length(area) != nrow(proj_X_ij) && length(area) != 1L) {
      stop("`area` should be of the same length as `nrow(newdata)` or of length 1.", call. = FALSE)
    }

    tmb_data$proj_X_threshold <- thresh$X_threshold
    tmb_data$area_i <- if (length(area) == 1L) rep(area, nrow(proj_X_ij)) else area
    tmb_data$proj_mesh <- proj_mesh
    tmb_data$proj_X_ij <- proj_X_ij
    tmb_data$proj_X_rw_ik <- proj_X_rw_ik
    tmb_data$proj_RE_indexes <- proj_RE_indexes
    tmb_data$proj_year <- make_year_i(nd[[object$time]])
    tmb_data$proj_lon <- newdata[[xy_cols[[1]]]]
    tmb_data$proj_lat <- newdata[[xy_cols[[2]]]]
    tmb_data$calc_se <- as.integer(se_fit)
    tmb_data$pop_pred <- as.integer(pop_pred)
    tmb_data$exclude_RE <- exclude_RE
    # tmb_data$calc_index_totals <- as.integer(!se_fit)
    # tmb_data$calc_cog <- as.integer(!se_fit)
    tmb_data$proj_spatial_index <- newdata$sdm_spatial_id
    tmb_data$proj_Zs <- sm$Zs
    tmb_data$proj_Xs <- sm$Xs
    if (!is.null(object$spatial_varying)) {
      if (!object$spatial_varying %in% names(newdata))
        stop("The `spatial_varying` column is missing from `newdata`.", call. = FALSE)
    }
    tmb_data$proj_z_i <- if (is.null(object$spatial_varying)) 0 else newdata[[object$spatial_varying]]

    epsilon_covariate <- rep(0, length(unique(newdata[[object$time]])))
    if (tmb_data$est_epsilon_model) {
      # covariate vector dimensioned by number of time steps
      time_steps <- unique(newdata[[object$time]])
      for (i in seq_along(time_steps)) {
        epsilon_covariate[i] <- newdata[newdata[[object$time]] == time_steps[i],
            object$epsilon_predictor, drop = TRUE][[1]]
      }
    }
    tmb_data$epsilon_predictor <- epsilon_covariate

    new_tmb_obj <- TMB::MakeADFun(
      data = tmb_data,
      parameters = get_pars(object),
      map = object$tmb_map,
      random = object$tmb_random,
      DLL = "sdmTMB",
      silent = TRUE
    )

    old_par <- object$model$par
    # need to initialize the new TMB object once:
    new_tmb_obj$fn(old_par)

    if (sims > 0) {
      if (!"jointPrecision" %in% names(object$sd_report) && !has_no_random_effects(object)) {
        message("Rerunning TMB::sdreport() with `getJointPrecision = TRUE`.")
        sd_report <- TMB::sdreport(object$tmb_obj, getJointPrecision = TRUE)
      } else {
        sd_report <- object$sd_report
      }
      if (has_no_random_effects(object)) {
        t_draws <- t(mvtnorm::rmvnorm(n = sims, mean = sd_report$par.fixed,
          sigma = sd_report$cov.fixed))
        row.names(t_draws) <- NULL
      } else {
        # if (!mvn_mle) {
          t_draws <- rmvnorm_prec(mu = new_tmb_obj$env$last.par.best,
            tmb_sd = sd_report, n_sims = sims)
        # } else {
        #   t_draws <- rmvnorm_prec_random(new_tmb_obj,
        #     tmb_sd = sd_report, n_sims = sims)
        # }
      }
      r <- apply(t_draws, 2L, new_tmb_obj$report)
    }
    if (!is.null(tmbstan_model)) {
      if (!"stanfit" %in% class(tmbstan_model))
        stop("tmbstan_model must be output from tmbstan::tmbstan().", call. = FALSE)
      t_draws <- extract_mcmc(tmbstan_model)
      r <- apply(t_draws, 2L, new_tmb_obj$report)
    }
    if (!is.null(tmbstan_model) || sims > 0) {
      .var <-  switch(sims_var,
        "est" = "proj_eta",
        "est_rf" = "proj_rf",
        "omega_s" = "proj_re_sp_st",
        "zeta_s" = "proj_re_sp_slopes",
        "epsilon_st" = "proj_re_st_vector",
        sims_var)
      out <- lapply(r, `[[`, .var)
      out <- do.call("cbind", out)
      rownames(out) <- nd[[object$time]] # for use in index calcs
      attr(out, "time") <- object$time
      return(out)
    }

    lp <- new_tmb_obj$env$last.par.best
    r <- new_tmb_obj$report(lp)

    if (isFALSE(pop_pred)) {
      nd$est <- r$proj_eta
      nd$est_non_rf <- r$proj_fe
      nd$est_rf <- r$proj_rf
      nd$omega_s <- r$proj_re_sp_st
      nd$zeta_s <- r$proj_re_sp_slopes
      nd$epsilon_st <- r$proj_re_st_vector
    }

    nd$sdm_spatial_id <- NULL
    nd$sdm_orig_id <- NULL

    obj <- new_tmb_obj

    if (se_fit) {
      sr <- TMB::sdreport(new_tmb_obj, bias.correct = FALSE)
      ssr <- summary(sr, "report")
      if (pop_pred) {
        proj_eta <- ssr[row.names(ssr) == "proj_fe", , drop = FALSE]
      } else {
        proj_eta <- ssr[row.names(ssr) == "proj_eta", , drop = FALSE]
      }
      row.names(proj_eta) <- NULL
      d <- as.data.frame(proj_eta)
      names(d) <- c("est", "se")
      nd$est <- d$est
      nd$est_se <- d$se
    }

    if ("sdmTMB_fake_year" %in% names(nd)) {
      nd <- nd[!nd$sdmTMB_fake_year,,drop=FALSE]
      nd$sdmTMB_fake_year <- NULL
    }
    if (fake_spatial_added) {
      for (i in 1:2) nd[[xy_cols[[i]]]] <- NULL
    }
    if (sdmTMB_fake_response) {
      nd[[response]] <- NULL
    }

  } else { # We are not dealing with new data:
    if (se_fit) {
      warning("Standard errors have not been implemented yet unless you ",
        "supply `newdata`. In the meantime you could supply your original data frame ",
        "to the `newdata` argument.", call. = FALSE)
    }
    nd <- object$data
    lp <- object$tmb_obj$env$last.par.best
    # object$tmb_obj$fn(lp) # call once to update internal structures?
    r <- object$tmb_obj$report(lp)

    nd$est <- r$eta_i
    # The following is not an error,
    # IID and RW effects are baked into fixed effects for `newdata` in above code:
    nd$est_non_rf <- r$eta_fixed_i + r$eta_rw_i + r$eta_iid_re_i
    nd$est_rf <- r$omega_s_A + r$epsilon_st_A_vec + r$zeta_s_A
    nd$omega_s <- r$omega_s_A
    nd$zeta_s <- r$zeta_s_A
    nd$epsilon_st <- r$epsilon_st_A_vec
    obj <- object
  }

  if (return_tmb_object)
    return(list(data = nd, report = r, obj = obj, fit_obj = object, pred_tmb_data = tmb_data))
  else {
    if (visreg_df) {
      # for visreg & related, return consistent objects with lm(), gam() etc.
      if (isTRUE(se_fit)) {
        return(list(fit = nd[,"est"], se.fit = nd[,"est_se"]))
      } else {
        return(as.numeric(nd[,"est"])) # return numeric vector
      }
    } else {
      # return dataframe by default
      return(nd)
    }
  }
}

# https://stackoverflow.com/questions/13217322/how-to-reliably-get-dependent-variable-name-from-formula-object
get_response <- function(formula) {
  tt <- terms(formula)
  vars <- as.character(attr(tt, "variables"))[-1] ## [1] is the list call
  response <- attr(tt, "response") # index of response var
  vars[response]
}

remove_9000 <- function(x) {
  as.package_version(paste0(
    strsplit(as.character(x), ".", fixed = TRUE)[[1]][1:3],
    collapse = "."
  ))
}

check_sdmTMB_version <- function(version) {
  if (remove_9000(utils::packageVersion("sdmTMB")) >
    remove_9000(version)) {
    warning(
      "The installed version of sdmTMB is newer than the version\n",
      "that was used to fit this model. It is possible new parameters\n",
      "have been added to the TMB model since you fit this model and\n",
      "that prediction will fail. We recommend you fit and predict\n",
      "from an sdmTMB model with the same version.",
      call. = FALSE
    )
  }
}

check_time_class <- function(object, newdata) {
  cls1 <- class(object$data[[object$time]])
  cls2 <- class(newdata[[object$time]])
  if (!identical(cls1, cls2)) {
    if (!identical(sort(c(cls1, cls2)), c("integer", "numeric"))) {
      stop(
        "Class of fitted time column (", cls1, ") does not match class of ",
        "`newdata` time column (", cls2 ,").",
        call. = FALSE)
    }
  }
}

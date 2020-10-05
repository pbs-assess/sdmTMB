# @importFrom stats predict
# @rdname predict
# @export

#' Predict from an sdmTMB model
#'
#' Can predict on the original data locations or on new data.
#'
#' @param object An object from [sdmTMB()].
#' @param newdata An optional new data frame. This should be a data frame with
#'   the same predictor columns as in the fitted data and a time column (if this
#'   is a spatiotemporal model) with the same name as in the fitted data. There
#'   should be predictor data for each year in the original data set.
#' @param se_fit Should standard errors on predictions at the new locations given by
#'   `newdata` be calculated? Warning: the current implementation can be slow for
#'   large data sets or high-resolution projections.
#' @param xy_cols A character vector of length 2 that gives the column names of
#'   the x and y coordinates in `newdata`.
#' @param return_tmb_object Logical. If `TRUE`, will include the TMB object in
#'   a list format output. Necessary for the [get_index()] or [get_cog()] functions.
#' @param area A vector of areas for survey grid cells. Only necessary if the
#'   output will be passed to [get_index()] or [get_cog()]. Should be the same length
#'   as the number of rows of `newdata`. Defaults to a sequence of 1s.
#' @param re_form `NULL` to specify individual-level predictions. `~0` or `NA`
#'   for population-level predictions. Note that unlike lme4 or glmmTMB, this
#'   only affects what the standard errors are calculated on if `se_fit = TRUE`.
#'   Otherwise, predictions at various levels are returned in all cases.
#' @param ... Not implemented.
#'
#' @return
#' If `return_tmb_object = FALSE`:
#' A data frame:
#' * `est`: Estimate in link space (everything is in link space)
#' * `est_non_rf`: Estimate from everything that isn't a random field
#' * `est_rf`: Estimate from all random fields combined
#' * `omega_s`: Spatial (intercept) random field that is constant through time
#' * `zeta_s`: Spatial slope random field
#' * `epsilon_st`: Spatiotemporal (intercept) random fields (could be
#'    independent draws each year or AR1)
#'
#' If `return_tmb_object = TRUE`:
#' A list:
#' * `data`: The data frame described above
#' * `report`: The TMB report on parameter values
#' * `obj`: The TMB object returned from the prediction run.
#' * `fit_obj`: The original TMB model object.
#'
#' You likely only need the `data` element as an end user. The other elements
#' are included for other functions.
#'
#' @export
#'
#' @examples
#' # We'll only use a small number of knots so this example runs quickly
#' # but you will likely want to use many more in applied situations.
#'
#' library(ggplot2)
#' d <- pcod
#' pcod_spde <- make_spde(d, c("X", "Y"), cutoff = 30) # a coarse mesh for example speed
#' m <- sdmTMB(
#'  data = d, formula = density ~ 0 + as.factor(year) + depth_scaled + depth_scaled2,
#'  time = "year", spde = pcod_spde, family = tweedie(link = "log")
#' )
#'
#' # Predictions at original data locations -------------------------------
#'
#' predictions <- predict(m)
#' head(predictions)
#'
#' predictions$resids <- residuals(m) # randomized quantile residuals
#' ggplot(predictions, aes(X, Y, col = resids)) + scale_colour_gradient2() +
#'   geom_point() + facet_wrap(~year)
#' hist(predictions$resids)
#' qqnorm(predictions$resids);abline(a = 0, b = 1)
#'
#' # Predictions onto new data --------------------------------------------
#'
#' predictions <- predict(m, newdata = qcs_grid)
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
#' \donttest{
#' # Visualizing a marginal effect ----------------------------------------
#' # Also demonstrates getting standard errors on population-level predictions
#'
#' nd <- data.frame(depth_scaled =
#'   seq(min(d$depth_scaled), max(d$depth_scaled), length.out = 100))
#' nd$depth_scaled2 <- nd$depth_scaled^2
#'
#' # You'll need at least one time element. If time isn't also a fixed effect
#' # then it doesn't matter what you pick:
#' nd$year <- 2003L
#' p <- predict(m, newdata = nd, se_fit = TRUE, re_form = NA)
#' ggplot(p, aes(depth_scaled, exp(est),
#'   ymin = exp(est - 1.96 * est_se), ymax = exp(est + 1.96 * est_se))) +
#'   geom_line() + geom_ribbon(alpha = 0.4)
#'
#' # Plotting marginal effect of a spline ---------------------------------
#'
#' m_gam <- sdmTMB(
#'  data = d, formula = density ~ 0 + as.factor(year) + s(depth_scaled, k = 3),
#'  time = "year", spde = pcod_spde, family = tweedie(link = "log")
#' )
#' nd <- data.frame(depth_scaled =
#'   seq(min(d$depth_scaled), max(d$depth_scaled), length.out = 100))
#' nd$year <- 2003L
#' p <- predict(m_gam, newdata = nd, se_fit = TRUE, re_form = NA)
#' ggplot(p, aes(depth_scaled, exp(est),
#'   ymin = exp(est - 1.96 * est_se), ymax = exp(est + 1.96 * est_se))) +
#'   geom_line() + geom_ribbon(alpha = 0.4)
#'
#' # Forecasting ----------------------------------------------------------
#'
#' # Forecasting using Eric Ward's hack with the `weights` argument.
#'
#' # Add on a fake year of data with the year to forecast:
#' dfake <- rbind(d, d[nrow(d), ])
#' dfake[nrow(dfake), "year"] <- max(d$year) + 1
#' tail(dfake) # last row is fake!
#'
#' weights <- rep(1, nrow(dfake))
#' weights[length(weights)] <- 0 # set last row weight to 0
#' dfake$year_factor <- dfake$year
#' dfake$year_factor[nrow(dfake)] <- max(d$year) # share fixed effect for last 2 years
#'
#' pcod_spde <- make_spde(dfake, c("X", "Y"), cutoff = 30)
#'
#' m <- sdmTMB(
#'   data = dfake, formula = density ~ 0 + as.factor(year_factor),
#'   ar1_fields = TRUE, # using an AR1 to have something to forecast with
#'   weights = weights,
#'   include_spatial = TRUE, # could also be `FALSE`
#'   time = "year", spde = pcod_spde, family = tweedie(link = "log")
#' )
#'
#' # Add a year to our grid:
#' qcs_grid$year_factor <- qcs_grid$year
#' grid2018 <- qcs_grid[qcs_grid$year == 2017L, ]
#' grid2018$year <- 2018L # `L` because `year` is an integer in the data
#' qcsgrid_forecast <- rbind(qcs_grid, grid2018)
#'
#' predictions <- predict(m, newdata = qcsgrid_forecast)
#' plot_map(predictions, "exp(est)") +
#'   scale_fill_viridis_c(trans = "log10")
#' plot_map(predictions, "epsilon_st") +
#'   scale_fill_gradient2()
#'
#' # Estimating local trends ----------------------------------------------
#'
#' pcod_spde <- make_spde(pcod, c("X", "Y"), cutoff = 25)
#' m <- sdmTMB(data = pcod, formula = density ~ depth_scaled + depth_scaled2,
#'   spde = pcod_spde, family = tweedie(link = "log"),
#'   spatial_trend = TRUE, time = "year", spatial_only = TRUE)
#' p <- predict(m, newdata = qcs_grid)
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

predict.sdmTMB <- function(object, newdata = NULL, se_fit = FALSE,
  xy_cols = c("X", "Y"), return_tmb_object = FALSE,
  area = 1, re_form = NULL, ...) {

  if (!missing(xy_cols)) {
    warning("argument `xy_cols` is deprecated; this information is already ",
    "in the output of the new `make_spde().", call. = FALSE)
  }
  if (!"xy_cols" %in% names(object$spde)) {
    warning("It looks like this model was fit with an older version of make_spde(). ",
    "Using `xy_cols`, but future versions of sdmTMB may not be compatible with this ",
    "model. Please update the make_spde() object and model fit.", call. = FALSE)
  } else {
    xy_cols <- object$spde$xy_cols
  }

  test <- suppressWarnings(tryCatch(object$tmb_obj$report(), error = function(e) NA))
  if (all(is.na(test))) object <- update_model(object)

  # from glmmTMB:
  pop_pred <- (!is.null(re_form) && ((re_form == ~0) || identical(re_form, NA)))

  tmb_data <- object$tmb_data
  tmb_data$do_predict <- 1L

  if (!is.null(newdata)) {
    if (any(!xy_cols %in% names(newdata)) && isFALSE(pop_pred))
      stop("`xy_cols` (the column names for the x and y coordinates) ",
        "are not in `newdata`. Did you miss specifying the argument ",
        "`xy_cols` to match your data?", call. = FALSE)

    if (object$time == "_sdmTMB_time") newdata[[object$time]] <- 0L
    if (!identical(class(object$data[[object$time]]), class(newdata[[object$time]])))
      stop("Class of fitted time column does not match class of `newdata` time column.",
        call. = FALSE)
    original_time <- sort(unique(object$data[[object$time]]))
    new_data_time <- sort(unique(newdata[[object$time]]))

    if (!all(new_data_time %in% original_time))
      stop("Some new time elements were found in `newdata`. ",
      "For now, make sure only time elements from the original dataset are present. If you would like to predict on new time elements, see the example hack with the `weights` argument in the help for `?predict.sdmTMB`.",
        call. = FALSE)

    if (!all(original_time %in% new_data_time)) {
      newdata[["sdmTMB_fake_year"]] <- FALSE
      missing_time_elements <- original_time[!original_time %in% new_data_time]
      nd2 <- do.call("rbind",
        replicate(length(missing_time_elements), newdata[1L,,drop=FALSE], simplify = FALSE))
      nd2[[object$time]] <- rep(missing_time_elements, each = 1L)
      nd2$sdmTMB_fake_year <- TRUE
      newdata <- rbind(newdata, nd2)
    }

    # If making population predictions (with standard errors), we don't need
    # to worry about space, so fill in dummy values if the user hasn't made any:
    fake_spatial_added <- FALSE
    if (pop_pred) {
      for (i in 1:2) {
        if (!xy_cols[[i]] %in% names(newdata)) {
          newdata[[xy_cols[[i]]]] <- mean(object$data[[xy_cols[[i]]]], na.rm = TRUE)
          fake_spatial_added <- TRUE
        }
      }
    }

    if (sum(is.na(new_data_time)) > 1)
      stop("There is at least one NA value in the time column. ",
        "Please remove it.", call. = FALSE)

    newdata$sdm_orig_id <- seq(1, nrow(newdata))
    fake_newdata <- unique(newdata[,xy_cols])
    fake_newdata[["sdm_spatial_id"]] <- seq(1, nrow(fake_newdata)) - 1L

    newdata <- base::merge(newdata, fake_newdata, by = xy_cols,
      all.x = TRUE, all.y = FALSE)
    newdata <- newdata[order(newdata$sdm_orig_id),, drop=FALSE]

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

    if (!"mgcv" %in% names(object)) object[["mgcv"]] <- FALSE
    proj_X_ij <- matrix(999)
    if (!object$mgcv) {
      proj_X_ij <- tryCatch({model.matrix(object$formula, data = nd)},
        error = function(e) NA)
    }
    if (object$mgcv || identical(proj_X_ij, NA)) {
      proj_X_ij <- mgcv::predict.gam(object$mgcv_mod, type = "lpmatrix", newdata = nd)
    }
    if (!is.null(object$time_varying))
      proj_X_rw_ik <- model.matrix(object$time_varying, data = nd)
    else
      proj_X_rw_ik <- matrix(0, ncol = 1, nrow = 1) # dummy

    tmb_data$proj_X_threshold <- thresh$X_threshold
    tmb_data$area_i <- if (length(area) == 1L && area[[1]] == 1) rep(1, nrow(proj_X_ij)) else area
    tmb_data$proj_mesh <- proj_mesh
    tmb_data$proj_X_ij <- proj_X_ij
    tmb_data$proj_X_rw_ik <- proj_X_rw_ik
    tmb_data$proj_year <- make_year_i(nd[[object$time]])
    tmb_data$proj_lon <- newdata[[xy_cols[[1]]]]
    tmb_data$proj_lat <- newdata[[xy_cols[[2]]]]
    tmb_data$calc_se <- as.integer(se_fit)
    tmb_data$pop_pred <- as.integer(pop_pred)
    tmb_data$calc_time_totals <- as.integer(!se_fit)
    tmb_data$proj_spatial_index <- newdata$sdm_spatial_id
    tmb_data$proj_t_i <- as.numeric(newdata[[object$time]])
    tmb_data$proj_t_i <- tmb_data$proj_t_i - mean(unique(tmb_data$proj_t_i)) # center on mean
    new_tmb_obj <- TMB::MakeADFun(
      data = tmb_data,
      parameters = object$tmb_obj$env$parList(),
      map = object$tmb_map,
      random = object$tmb_random,
      DLL = "sdmTMB",
      silent = TRUE
    )

    old_par <- object$model$par
    # need to initialize the new TMB object once:
    new_tmb_obj$fn(old_par)
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
    lp <- object$tmb_obj$env$last.par
    # object$tmb_obj$fn(lp) # call once to update internal structures?
    r <- object$tmb_obj$report(lp)

    nd$est <- r$eta_i
    # Following is not an error: rw effects baked into fixed effects for new data in above code:
    nd$est_non_rf <- r$eta_fixed_i + r$eta_rw_i
    nd$est_rf <- r$omega_s_A + r$epsilon_st_A_vec + r$omega_s_trend_A
    nd$omega_s <- r$omega_s_A
    nd$zeta_s <- r$omega_s_trend_A
    nd$epsilon_st <- r$epsilon_st_A_vec
    obj <- object
  }

  if (return_tmb_object)
    return(list(data = nd, report = r, obj = obj, fit_obj = object))
  else
    return(nd)
}

# https://stackoverflow.com/questions/13217322/how-to-reliably-get-dependent-variable-name-from-formula-object
get_response <- function(formula) {
  tt <- terms(formula)
  vars <- as.character(attr(tt, "variables"))[-1] ## [1] is the list call
  response <- attr(tt, "response") # index of response var
  vars[response]
}

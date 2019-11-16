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
#'
#' pcod_spde <- make_spde(d$X, d$Y, n_knots = 50) # just 50 for example speed
#' m <- sdmTMB(
#'  data = d, formula = density ~ 0 + as.factor(year) + depth_scaled + depth_scaled2,
#'  time = "year", spde = pcod_spde, family = tweedie(link = "log")
#' )
#'
#' # Predictions at original data locations:
#' predictions <- predict(m)
#' head(predictions)
#'
#' predictions$resids <- residuals(m) # randomized quantile residuals
#' ggplot(predictions, aes(X, Y, col = resids)) + scale_colour_gradient2() +
#'   geom_point() + facet_wrap(~year)
#' hist(predictions$resids)
#' qqnorm(predictions$resids);abline(a = 0, b = 1)
#'
#' # Predictions onto new data:
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
#'   scale_fill_viridis_c(trans = "sqrt")
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
#' # Spatial trend example:
#' pcod_spde <- make_spde(d$X, d$Y, n_knots = 100)
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
#'
#' # Example with standard errors on new location predictions.
#' # Note this example models presence/absence.
#' pcod_2017 <- d[d$year == 2017, ]
#' pcod_spde <- make_spde(pcod_2017$X, pcod_2017$Y, n_knots = 75)
#' m2017 <- sdmTMB(
#'   pcod_2017, present ~ 0 + depth_scaled + depth_scaled2,
#'   time = "year", spde = pcod_spde, family = binomial(link = "logit"),
#'   silent = FALSE
#' )
#' }

predict.sdmTMB <- function(object, newdata = NULL, se_fit = FALSE,
  xy_cols = c("X", "Y"), return_tmb_object = FALSE,
  area = 1, ...) {

  tmb_data <- object$tmb_data
  tmb_data$do_predict <- 1L
  TMB::openmp(1L)

  if (!is.null(newdata)) {
    if (any(!xy_cols %in% names(newdata)))
      stop("`xy_cols` (the column names for the x and y coordinates) ",
        "are not in `newdata`. Did you miss specifying the argument ",
        "to match your data?", call. = FALSE)

    if (object$time == "_sdmTMB_time") newdata[[object$time]] <- 0L
    original_time <- sort(unique(object$data[[object$time]]))
    new_data_time <- sort(unique(newdata[[object$time]]))
    if (!identical(original_time, new_data_time))
      stop("For now, all of the time elements in the original data set must ",
      "be identical to the time elements in the `newdata` data set ",
      "but they are not.", call. = FALSE)
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

    nd <- newdata
    response <- get_response(object$formula)
    if (!response %in% names(nd)) nd[[response]] <- 0 # fake for model.matrix
    proj_X_ij <- model.matrix(object$formula, data = nd)
    if (!is.null(object$time_varying))
      proj_X_rw_ik <- model.matrix(object$time_varying, data = nd)
    else
      proj_X_rw_ik <- matrix(0, ncol = 1, nrow = 1) # dummy


    # mf   <- model.frame(object$formula, data = nd)
    # offset <- as.vector(model.offset(mf))
    # if (is.null(offset)) offset <- rep(0, nrow(nd))

    # tmb_data$offset_i <- offset
    tmb_data$area_i <- if (length(area) == 1L && area[[1]] == 1) rep(1, nrow(proj_X_ij)) else area
    tmb_data$proj_mesh <- proj_mesh
    tmb_data$proj_X_ij <- proj_X_ij
    tmb_data$proj_X_rw_ik <- proj_X_rw_ik
    tmb_data$proj_year <- make_year_i(nd[[object$time]])
    tmb_data$proj_lon <- newdata[[xy_cols[[1]]]]
    tmb_data$proj_lat <- newdata[[xy_cols[[2]]]]
    tmb_data$calc_se <- as.integer(se_fit)
    tmb_data$calc_time_totals <- 1L # for now (always on)
    tmb_data$proj_spatial_index <- newdata$sdm_spatial_id
    tmb_data$proj_t_i <- as.numeric(newdata[[object$time]])
    tmb_data$proj_t_i <- tmb_data$proj_t_i - mean(unique(tmb_data$proj_t_i)) # center on mean
    new_tmb_obj <- TMB::MakeADFun(
      data = tmb_data,
      parameters = object$tmb_params,
      map = object$tmb_map,
      random = object$tmb_random,
      DLL = "sdmTMB",
      silent = TRUE
    )

    old_par <- object$model$par
    # need to initialize the new TMB object once:
    new_tmb_obj$fn(old_par)
    lp <- new_tmb_obj$env$last.par

    r <- new_tmb_obj$report(lp)

    nd$est <- r$proj_eta
    nd$est_non_rf <- r$proj_fe
    nd$est_rf <- r$proj_rf
    nd$omega_s <- r$proj_re_sp_st
    nd$zeta_s <- r$proj_re_sp_slopes
    nd$epsilon_st <- r$proj_re_st_vector

    nd$sdm_spatial_id <- NULL
    nd$sdm_orig_id <- NULL

    obj <- new_tmb_obj

    if (se_fit) {
      sr <- TMB::sdreport(new_tmb_obj, bias.correct = FALSE)
      ssr <- summary(sr, "report")
      proj_eta <- ssr[row.names(ssr) == "proj_eta", , drop = FALSE]
      row.names(proj_eta) <- NULL
      d <- as.data.frame(proj_eta)
      names(d) <- c("est", "se")
      nd$est_se <- d$se
    }
  } else {
    if (se_fit) {
      warning("Standard errors have not been implemented yet unless you ",
        "supply `newdata`. In the meantime you could supply your original data frame ",
        "to the `newdata` argument.", call. = FALSE)
    }
    nd <- object$data
    lp <- object$tmb_obj$env$last.par
    r <- object$tmb_obj$report(lp)

    nd$est <- r$eta_i
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

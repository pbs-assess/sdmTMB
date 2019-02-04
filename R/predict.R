# @importFrom stats predict
# @rdname predict
# @export

#' Predict from an sdmTMB model
#'
#' Can predict on the original data locations or onto new data.
#'
#' @param object An object from [sdmTMB()].
#' @param newdata An optional new data frame. This should be a single set of
#'   spatial locations. These locations will be expanded to cover all the years
#'   in the original data set. Eventually `newdata` will be more flexible.
#' @param se_fit Should standard errors on predictions at the new locations given by
#'   `newdata` be calculated? Warning: the current implementation can be slow for
#'   large data sets or high-resolution projections.
#' @param xy_cols A character vector of length 2 that gives the column names of
#'   the x and y coordinates in `newdata`.
#' @param ... Not implemented.
#'
#' @return
#' A data frame. TODO details.
#' @export
#'
#' @examples
#' # We'll only use a small number of knots so this example runs quickly
#' # but you will likely want to use many more (depending on your data).
#'
#' pcod_spde <- make_spde(pcod$X, pcod$Y, n_knots = 100)
#' m <- sdmTMB(
#'  pcod, density ~ 0 + as.factor(year) + depth_scaled + depth_scaled2,
#'  time = "year", spde = pcod_spde, family = tweedie(link = "log"),
#'  silent = FALSE
#' )
#'
#' # Predictions at original data locations:
#' predictions <- predict(m)$data
#' cols <- c("year", "X", "Y", "est", "est_fe",
#'   "est_re_s", "est_re_st", "s_i")
#' head(predictions[,cols])
#'
#' predictions$resids <- residuals(m) # randomized quantile residuals
#' library(ggplot2)
#' ggplot(predictions, aes(X, Y, col = resids)) + scale_colour_gradient2() +
#'   geom_point() + facet_wrap(~year)
#' hist(predictions$resids)
#' qqnorm(predictions$resids);abline(a = 0, b = 1)
#'
#' # Predictions onto new data:
#' predictions <- predict(m, newdata = qcs_grid)$data
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
#' plot_map(predictions, "exp(est_fe)") +
#'   ggtitle("Prediction (fixed effects only)") +
#'   scale_fill_viridis_c(trans = "sqrt")
#'
#' plot_map(predictions, "est_re_s") +
#'   ggtitle("Spatial random effects only") +
#'   scale_fill_gradient2()
#'
#' plot_map(predictions, "est_re_st") +
#'   ggtitle("Spatiotemporal random effects only") +
#'   scale_fill_gradient2()
#'
#' \donttest{
#' # Example with standard errors on new location predictions.
#' # Note this example models presence/absence.
#' pcod_2017 <- pcod[pcod$year == 2017, ]
#' pcod_spde <- make_spde(pcod_2017$X, pcod_2017$Y, n_knots = 75)
#' m2017 <- sdmTMB(
#'   pcod_2017, present ~ 0 + depth_scaled + depth_scaled2,
#'   time = "year", spde = pcod_spde, family = binomial(link = "logit"),
#'   silent = FALSE
#' )
#'
#' # Predictions at new data locations with standard errors.
#' # Note that this can currently be quite slow on large data sets.
#' predictions <- predict(m2017, newdata = qcs_grid, se_fit = TRUE)$data
#'
#' plot_map(predictions, "plogis(est)") +
#'   scale_fill_gradient2(midpoint = 0.5) +
#'   ggtitle("Predictions")
#'
#' plot_map(predictions, "est_se") +
#'   scale_fill_viridis_c() +
#'   ggtitle("Prediction standard error")
#'
#' plot_map(predictions, "plogis(est + 2 * est_se)") +
#'   scale_fill_gradient2(midpoint = 0.5) +
#'   ggtitle("Prediction upper 95% CI")
#'
#' plot_map(predictions, "plogis(est - 2 * est_se)") +
#'   scale_fill_gradient2(midpoint = 0.5) +
#'   ggtitle("Prediction lower 95% CI")
#' }

predict.sdmTMB <- function(object, newdata = NULL, se_fit = FALSE,
  xy_cols = c("X", "Y"), ...) {

  tmb_data <- object$tmb_data
  tmb_data$do_predict <- 1L

  if (!is.null(newdata)) {
    proj_mesh <- INLA::inla.spde.make.A(object$spde$mesh,
      loc = as.matrix(newdata[, xy_cols]))

    # expand for time units:
    original_time <- sort(unique(object$data[[object$time]]))
    nd <- do.call("rbind",
      replicate(length(original_time), newdata, simplify = FALSE))
    nd[[object$time]] <- rep(original_time, each = nrow(newdata))

    nd[[get_response(object$formula)]] <- 0 # fake for model.matrix
    proj_X_ij <- model.matrix(object$formula, data = nd)

    tmb_data$proj_mesh <- proj_mesh
    tmb_data$proj_X_ij <- proj_X_ij
    tmb_data$proj_year <- as.integer(as.factor(as.character(nd$year))) - 1L
    tmb_data$calc_se <- as.integer(se_fit)
    tmb_data$calc_time_totals <- 1L # for now

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
    # new_tmb_obj$env$parList()

    r <- new_tmb_obj$report(lp)

    nd$est <- r$proj_eta
    nd$est_fe <- r$proj_fe
    nd$est_re_s <- r$proj_re_sp
    nd$est_re_st <- r$proj_re_st_vector
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
    nd$est_fe <- r$eta_fixed_i

    # Spatial REs:
    get_omegas <- function(x) r$omega_s[x]
    nd$est_re_s <- mapply(get_omegas, x = object$tmb_data$s_i + 1)

    # Sp-temp REs:
    epsilon_st_mat <- object$tmb_obj$report(lp)$epsilon_st
    get_eps <- function(x, y) epsilon_st_mat[x, y]
    eps <- mapply(get_eps,
      x = object$tmb_data$s_i + 1,
      y = object$tmb_data$year_i + 1)
    nd$est_re_st <- eps

    nd$s_i <- object$tmb_data$s_i + 1
    obj <- object
  }

  list(data = nd, report = r, obj = obj)
}

# https://stackoverflow.com/questions/13217322/how-to-reliably-get-dependent-variable-name-from-formula-object
get_response <- function(formula) {
  tt <- terms(formula)
  vars <- as.character(attr(tt, "variables"))[-1] ## [1] is the list call
  response <- attr(tt, "response") # index of response var
  vars[response]
}

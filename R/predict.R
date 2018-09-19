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
#' pcod_spde <- make_spde(pcod$X, pcod$Y, n_knots = 30)
#' m <- sdmTMB(
#'  pcod, density ~ 0 + as.factor(year) + depth_scaled + depth_scaled2,
#'  time = "year", spde = pcod_spde, family = tweedie(link = "log"),
#'  silent = TRUE
#' )
#' predictions <- predict(m)
#' cols <- c("year", "X", "Y", "prediction", "prediction_fe",
#'   "prediction_re", "s_i")
#' head(predictions[,cols])
#'
#' predictions <- predict(m, newdata = qcs_grid)
#' library(ggplot2)
#' plot_map <- function(dat, column = "prediction") {
#'   ggplot(dat, aes_string("X", "Y", fill = column)) +
#'     geom_raster() +
#'     facet_wrap(~year) +
#'     coord_fixed()
#' }
#' plot_map(predictions, "exp(prediction)") +
#'   scale_fill_viridis_c(trans = "sqrt") +
#'   ggtitle("Prediction (fixed effects + random effects)")
#' plot_map(predictions, "exp(prediction_fe)") +
#'   ggtitle("Prediction (fixed effects only)") +
#'   scale_fill_viridis_c(trans = "sqrt")
#' plot_map(predictions, "prediction_re") +
#'   ggtitle("Random effects only") +
#'   scale_fill_gradient2()

predict.sdmTMB <- function(object, newdata = NULL, xy_cols = c("X", "Y"), ...) {

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

    new_tmb_obj <- TMB::MakeADFun(
      data = tmb_data,
      parameters = object$tmb_params,
      random = object$random,
      DLL = "sdmTMB",
      silent = TRUE
    )

    old_par <- object$model$par
    # need to initialize the new TMB object once:
    new_tmb_obj$fn(old_par)
    lp <- new_tmb_obj$env$last.par

    nd$prediction <- new_tmb_obj$report(lp)$proj_eta
    nd$prediction_fe <- new_tmb_obj$report(lp)$proj_fe
    nd$prediction_re <- new_tmb_obj$report(lp)$proj_re_st_vector
  } else {
    nd <- object$data
    lp <- object$tmb_obj$env$last.par
    nd$prediction <- object$tmb_obj$report(lp)$eta_i
    nd$prediction_fe <- object$tmb_obj$report(lp)$linear_predictor_i
    epsilon_st_mat <- object$tmb_obj$report(lp)$epsilon_st
    get_eps <- function(x, y) epsilon_st_mat[x, y]
    eps <- mapply(get_eps,
      x = object$tmb_data$s_i + 1,
      y = object$tmb_data$year_i + 1)
    nd$prediction_re <- eps
    nd$s_i <- object$tmb_data$s_i + 1
  }
  nd
}

# https://stackoverflow.com/questions/13217322/how-to-reliably-get-dependent-variable-name-from-formula-object
get_response <- function(formula) {
  tt <- terms(formula)
  vars <- as.character(attr(tt, "variables"))[-1] ## [1] is the list call
  response <- attr(tt, "response") # index of response var
  vars[response]
}


ll_gaussian <- function(object, withheld_y, withheld_est) {
  dispersion <- exp(object$model$par[["ln_phi"]])
  dnorm(x = withheld_y, mean = withheld_est, sd = dispersion, log = TRUE)
}

ll_sdmTMB <- function(object, withheld_y, withheld_est, ...) {
  family_func <- switch(object$family$family,
    gaussian = ll_gaussian
    # binomial = nll_binomial,
    # tweedie  = nll_tweedie
  )
  family_func(object, withheld_y, withheld_est)
}

#' Save log likelihoods of k-fold cross-validation for sdmTMB models
#'
#' @param all_data A data frame.
#' @param formula Model formula.
#' @param family The family and link. Currently supports [gaussian()].
#' @param time Name of the time column.
#' @param x_coord Name of the column with X coordinate (as numeric vector).
#' @param y_coord Name of the column with Y coordinate (as numeric vector).
#' @param k_folds Number of folds.
#' @param fold_ids Optional input name of column containing user chosen fold ids.
#' @param n_knots The number of knots.
#' @param plot_spde Logical for whether or not to produce a plot of each mesh.
#' @param ... All other arguments required to run sdmTMB model with the exception of:
#'            [data] and [spde] which are redefined for each fold within the function.
#'
#' @return
#' @export
#'
#' @examples
#' d <- subset(pcod, year >= 2011) # subset for example speed
#'
#' # Gaussian
#' pcod_gaus <- subset(d, density > 0)
#' kfold <- loglik_cv(pcod_gaus, formula = log(density) ~ 0 + as.factor(year) + depth_scaled + depth_scaled2,
#'                  time = "year", x_coord = "X", y_coord = "Y",
#'                  n_knots = 30, k_folds = 3)
#' sum(kfold$data$cv_loglik)

loglik_cv <- function(all_data, formula, family, time = "year", x_coord = "X", y_coord = "Y", k_folds = 10, fold_ids = NULL, n_knots = NULL, plot_spde = TRUE, ...) {

  all_data <- as.data.frame(all_data)
  all_data$X <- all_data[[x_coord]]
  all_data$Y <- all_data[[y_coord]]

  # split data by 'time' from sdmTMB model arguments
  split_time <- all_data[[time]]

  # add column of fold_ids stratified across time steps
  if (is.null(fold_ids)) {
   dd <- lapply(split(all_data, split_time), function(x) {
     obs <- nrow(x)
     i <- obs / k_folds
     i <- round(c(0, i * seq(1, (k_folds - 1)), obs))
     times <- i[-1] - i[-length(i)]
     group <- c()
     for (j in 1:(length(times))) {
       group <- c(group, rep(j, times = times[j]))
     }
     r <- order(runif(obs))
     x$fold_ids <- group[r]
     x
   })
   d <- do.call(rbind, dd)
  } else {
   d <- all_data
   d$fold_ids <- all_data[[fold_ids]]
  }

  # model data k times for for k-1 folds
  out <- lapply(seq_len(k_folds), function(k) {
    d_fit <- d[d$fold_ids != k, , drop = FALSE]
    d_withheld <- d[d$fold_ids == k, , drop = FALSE]

    # build mesh for training data
    # FIXME: should we let the user set more parameters for the inla mesh?
    d_fit_spde <- make_spde(d_fit$X, d_fit$Y, n_knots = n_knots)

    if (plot_spde)
      plot_spde(d_fit_spde)

    # run model
    object <- sdmTMB(data = d_fit, formula = formula, time = time, spde = d_fit_spde, ...)

    # FIXME: Error in match.call(definition = def, call = def.call) :
    #         ... used in a situation where it does not exist
    # this is a line in get_args(), currently disabled
    # Chris working on repair to gfutilities::get.args()

    # predict for withheld data
    predicted <- predict(object, newdata = d_withheld)
    cv_data <- d_withheld
    cv_data$cv_est <- predicted$data$est

    response <- get_response(object$formula)
    withheld_y <- predicted$data[[response]]
    withheld_est <- predicted$data$est

    # calculate log likelihood for each withheld observation
    cv_data$cv_loglik <- ll_sdmTMB(object, withheld_y, withheld_est)

    list(data = cv_data, model = object)
   })
  data <- lapply (out, function(x) x$data)
  data <- do.call(rbind, data)
  models <- lapply (out, function(x) x$model)
  list(data = data, models = models)
}


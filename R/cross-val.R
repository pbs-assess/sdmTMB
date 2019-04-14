ll_gaussian <- function(object, withheld_y, withheld_mu) {
  .sd <- exp(object$model$par[["ln_phi"]])
  stats::dnorm(x = withheld_y, mean = withheld_mu, sd = .sd, log = TRUE)
}

ll_tweedie <- function(object, withheld_y, withheld_mu) {
  p <- stats::plogis(object$model$par[["thetaf"]]) + 1
  phi <- exp(object$model$par[["ln_phi"]])
  log(tweedie::dtweedie(y = withheld_y, power = p, mu = withheld_mu, phi = phi))
}

ll_sdmTMB <- function(object, withheld_y, withheld_mu) {
  family_func <- switch(object$family$family,
    gaussian = ll_gaussian,
    tweedie = ll_tweedie
  )
  family_func(object, withheld_y, withheld_mu)
}

#' Save log likelihoods of k-fold cross-validation for sdmTMB models
#'
#' @param formula Model formula.
#' @param data A data frame.
#' @param time Name of the time column.
#' @param x Name of the column with X coordinates.
#' @param y Name of the column with Y coordinates.
#' @param k_folds Number of folds.
#' @param fold_ids Optional input name of column containing user chosen fold
#'   ids.
#' @param n_knots The number of knots.
#' @param spde_function A function that takes 3 arguments (x, y, n_knots) and
#'   returns a list structure that matches output of [make_spde()].
#' @param plot_spde Logical for whether or not to produce a plot of each mesh.
#' @param ... All other arguments required to run sdmTMB model with the
#'   exception of `data` and `spde` which are redefined for each fold within the
#'   function.
#'
#' @export
#'
#' @examples
#' d <- subset(pcod, year >= 2011) # subset for example speed
#' x <- sdmTMB_cv(
#'   formula = density ~ 0 + depth_scaled + depth_scaled2 + as.factor(year),
#'   d, family = tweedie(link = "log"), time = "year", x = "X", y = "Y",
#'   n_knots = 30, k_folds = 3
#' )
#' x$sum_loglik
#' str(x$data)
#' x$models[[1]]
#' x$models[[2]]
#' \donttest{
#' # Parallel:
#' library(future)
#' plan(multiprocess)
#' x <- sdmTMB_cv(
#'   formula = density ~ 0 + depth_scaled + depth_scaled2 + as.factor(year),
#'   d, family = tweedie(link = "log"), time = "year", x = "X", y = "Y",
#'   n_knots = 30, k_folds = 4
#' )
#' }
sdmTMB_cv <- function(formula, data, time = "year", x = "X", y = "Y",
                      k_folds = 10, fold_ids = NULL, n_knots = NULL,
                      spde_function = make_spde, plot_spde = FALSE, ...) {
  data[["_sdm_order_"]] <- seq_len(nrow(data))

  # add column of fold_ids stratified across time steps
  if (is.null(fold_ids)) {
    dd <- lapply(split(data, data[[time]]), function(x) {
      obs <- nrow(x)
      i <- obs / k_folds
      i <- round(c(0, i * seq(1, (k_folds - 1)), obs))
      times <- i[-1] - i[-length(i)]
      group <- c()
      for (j in seq_along(times))
        group <- c(group, rep(j, times = times[j]))
      r <- order(stats::runif(obs))
      x$cv_fold <- group[r]
      x
    })
    data <- do.call(rbind, dd)
    fold_ids <- "cv_fold"
  }

  # model data k times for k-1 folds
  out <- future.apply::future_lapply(seq_len(k_folds), function(k) {
  # out <- lapply(seq_len(k_folds), function(k) {
    d_fit <- data[data[[fold_ids]] != k, , drop = FALSE]
    d_withheld <- data[data[[fold_ids]] == k, , drop = FALSE]

    # build mesh for training data
    d_fit_spde <- spde_function(d_fit[[x]], d_fit[[y]], n_knots = n_knots)

    # run model
    object <- sdmTMB(data = d_fit, formula = formula, time = time, spde = d_fit_spde, ...)

    # predict for withheld data
    predicted <- predict(object, newdata = d_withheld, xy_cols = c(x, y))
    cv_data <- d_withheld
    cv_data$cv_predicted <- object$family$linkinv(predicted$data$est)
    response <- get_response(object$formula)
    withheld_y <- predicted$data[[response]]
    withheld_mu <- cv_data$cv_predicted

    # calculate log likelihood for each withheld observationn
    cv_data$cv_loglik <- ll_sdmTMB(object, withheld_y, withheld_mu)
    list(data = cv_data, object = object)
  })
  models <- lapply(out, `[[`, "object")
  data <- lapply(out, `[[`, "data")
  data <- do.call(rbind, data)
  data <- data[order(data[["_sdm_order_"]]), , drop = FALSE]
  data[["_sdm_order_"]] <- NULL
  row.names(data) <- NULL

  if (plot_spde) {
    op <- graphics::par(
      mfrow = c(ceiling(sqrt(k_folds)), ceiling(sqrt(k_folds))),
      mar = c(1, 1, 1, 1)
    )
    for (i in seq_len(k_folds)){
      plot_spde(models[[i]]$spde)
    }
    graphics::par(op)
  }

  list(data = data, models = models, sum_loglik = sum(data$cv_loglik))
}

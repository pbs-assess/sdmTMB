# NOTE: dataframe must have coordinate data labeled as capital X and capital Y

ll_gaussian <- function(object, residuals) {
  dispersion <- exp(object$model$par[["ln_phi"]])
  # FIXME: This shouldn't be the residuals. The testing data
  # should be x, and the prediction should be mean.
  dnorm(x = residuals, mean = mean(residuals), sd = dispersion, log = TRUE)
}

ll_sdmTMB <- function(object, residuals, ...) {
  family_func <- switch(object$family$family,
    gaussian = ll_gaussian
    # binomial = nll_binomial,
    # tweedie  = nll_tweedie
  )
  family_func(object, residuals)
}

# FIXME: Make sure to add Roxygen documentation above for this function.
loglik_cv <- function(all_data, time = "year", k_folds = 10, fold_id = NULL, n_knots = NULL, ...) {

  # split data by 'time' from sdmTMB model arguments
  split_time <- all_data[[time]]

  # add column of fold_ids stratified across time steps
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
  # FIXME: Do we need the as.data.frame here?
  # If we do then we should adjust the code above so that we do not need it.
  d <- as.data.frame(do.call(rbind, dd))

  # model data k times for for k-1 folds
  out <- lapply(seq_len(k_folds), function(k) {
    d_fit <- d[d$fold_ids != k, , drop = FALSE]
    d_withheld <- d[d$fold_ids == k, , drop = FALSE]

    # build mesh for training data
    # FIXME: we should have the user pass the x and y column names so this can work with
    # any column names. Also, we should let the user pass a custom inla mesh.
    d_fit_spde <- make_spde(d_fit$X, d_fit$Y, n_knots = n_knots)

    # run model
    object <- sdmTMB(data = d_fit, spde = d_fit_spde, time = time, ...)

    # FIXME: Error in match.call(definition = def, call = def.call) :
    #         ... used in a situation where it does not exist
    # this is a line in get_args(), currently disabled
    # Chris working on repair to gfutilities::get.args()

    # predict for withheld data
    predicted <- predict(object, newdata = d_withheld)
    cv_data <- d_withheld
    cv_data$cv_est <- predicted$data$est

    # FIXME: We don't need to calculate the residuals.
    # calculate residuals
    response <- object$formula[[2]]
    residuals <- predicted$data[[response]] - predicted$data$est

    # calculate log likelihood for each withheld observation
    cv_data$cv_loglik <- ll_sdmTMB(object, residuals)
    cv_data

    list(data = cv_data, model = object)
  })
  data <- lapply (out, function(x) x$data)
  data <- do.call(rbind, data)
  models <- lapply (out, function(x) x$model)
  list(data = data, models = models)
}

# d_trawl <- readRDS("~/github/dfo/gfranges/analysis/tmb-sensor-explore/sensor-data-processed.rds")
# out <- loglik_cv(d_trawl, n_knots = 30, k_folds = 3, formula = temperature_c ~ 0 + as.factor(year))
# sum(out$data$nll)

# pcod <- load("~/github/dfo/sdmTMB/data/pcod.rda")

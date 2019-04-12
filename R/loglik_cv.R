# NOTE: dataframe must have coordinate data labeled as capital X and capital Y

nll_gaussian <- function(object, residuals) {
  dispersion <- exp(object$model$par[["ln_phi"]])
  nll <- -(dnorm(x = residuals, mean = mean(residuals), sd = dispersion, log = TRUE))
  nll
}

nll_sdmTMB <- function(object, residuals, ...) {
  family_func <- switch(object$family$family,
    gaussian = nll_gaussian
    # binomial = nll_binomial,
    # tweedie  = nll_tweedie
  )
  family_func(object, residuals)
}

loglik_cv <- function(all_data, time = "year", k_folds = 10, fold_id = NULL, n_knots = NULL, ...) {

  # split data by 'time' from sdmTMB model arguments
  split_time <- all_data[[time]]
  k <- k_folds

  # add column of fold_ids stratified across time steps
  dd <- lapply(split(all_data, split_time), function(x) {
    obs <- nrow(x)
    i <- obs / k
    i <- round(c(0, i * 1:(k - 1), obs))
    times <- i[-1] - i[-length(i)]
    group <- c()
    for (j in 1:(length(times))) {
      group <- c(group, rep(j, times = times[j]))
    }
    r <- order(runif(obs))
    x$fold_ids <- group[r]
    x
  })
  d <- as.data.frame(do.call(rbind, dd))

  # model data k times for for k-1 folds
  out <- lapply(1:k_folds, function(k) {
    d_fit <- d[d$fold_ids != k, , drop = FALSE]
    d_withheld <- d[d$fold_ids == k, , drop = FALSE]

    # build mesh for training data
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

    # calculate residuals
    response <- object$formula[[2]]
    residuals <- predicted$data[[response]] - predicted$data$est

    # calculate negative loglikelihood for each withheld observation
    cv_data$nll <- nll_sdmTMB(object, residuals)
    cv_data

    # would be nice to save model object too, but not sure how to
    # list(cv_data, object)
  })
  data <- do.call(rbind, out)
}

# d_trawl <- readRDS("~/github/dfo/gfranges/analysis/tmb-sensor-explore/sensor-data-processed.rds")
# out <- loglik_cv(d_trawl, n_knots = 30, k_folds = 3, formula = temperature_c ~ 0 + as.factor(year))
# sum(out$nll)

# pcod <- load("~/github/dfo/sdmTMB/data/pcod.rda")

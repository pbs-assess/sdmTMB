# pcod <- load("~/github/dfo/sdmTMB/data/pcod.rda")
d_trawl <- readRDS("~/github/dfo/gfranges/analysis/tmb-sensor-explore/sensor-data-processed.rds")
# NOTE: input data must have coordinate data labeled as capital X and capital Y


loglik_cv <- function(all_data, time = "year", k_folds = 10, fold_id = NULL, n_knots = NULL, ...) {
  #browser()
  # split data by 'time' from sdmTMB model arguments
  split_time <- all_data[[time]]
  k <- k_folds

  # add column of fold_ids stratified across time steps
  dd <- lapply(split(all_data, split_time), function(x) {
           obs <- nrow(x)
           i <- obs/k
           i <- round(c(0, i * 1:(k - 1), obs))
           times <- i[-1] - i[-length(i)]
           group <- c()
           for (j in 1:(length(times))) {
             group <- c(group, rep(j, times = times[j]))
           }
           r <- order(runif(obs))
           x$fold_ids <- group[r]
           x
          }
        )
  d <- as.data.frame(do.call(rbind, dd))

  # model data from k-1 folds
  for (i in seq_len(k_folds)) {
    d_fit <-d[d$fold_ids != i, , drop =FALSE]
    d_withheld <- d[d$fold_ids == i, , drop =FALSE]

    # build mesh for training data
    if (is.null(n_knots)) {
      n_time <- length(unique(split_time))
      auto_knots <- nrow(data)/n_time/2
      d_fit_spde <- sdmTMB::make_spde(d_fit$X, d_fit$Y, n_knots = auto_knots)
      } else {
      d_fit_spde <- sdmTMB::make_spde(d_fit$X, d_fit$Y, n_knots = n_knots)
      }
    browser()
    # run model
    m_fold <- sdmTMB::sdmTMB(data = d_fit, spde = d_fit_spde, time = time, ...)

#FIXME: Error in match.call(definition = def, call = def.call) :
#  ... used in a situation where it does not exist

    # predict for withheld data
    predicted <- predict(m_fold, newdata = d_withheld)

  }
}


loglik_cv(k_folds = 10, n_knots = 100, all_data = d_trawl, time = "year",
  formula = temperature_c ~ 0 + as.factor(year),
  time_varying = ~ 0 + depth_scaled,
  ar1_fields = FALSE,
  include_spatial = FALSE,
  family = gaussian(link = "identity"))

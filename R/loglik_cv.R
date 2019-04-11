data <- pcod
time = "year"

loglik_cv <- function(k_folds = 10, fold_id = NULL, ...) {

  # split data by 'time' from sdmTMB model arguments
  split_time <- data[[time]]

  # add column of fold_ids stratified across time steps
  dd <- lapply(split(data, split_time), function(x) {
           obs <- nrow(x)
           k <- k_folds
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
  do.call(rbind, dd)


  m_fold <- sdmTMB(     )


}

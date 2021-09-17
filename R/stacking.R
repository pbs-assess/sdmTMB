#' Perform stacking with log scores on sdmTMB_cv objects
#'
#' This approach is described in Yao et al. 2018, Bayesian Analysis. The general
#' method minimizes (or maximizes) some quantity across models. For simple models
#' with normal error, this may be the RMSE, but other approaches include log scores.
#' We adopt the latter here, where log scores are used to generate the stacking of
#' predictive distributions
#'
#' @param model_list A collection of models fit with sdmTMB_cv, to
#' generate estimates of predictive densities
#' @param include_folds An optional vector specifiying which folds to include
#' in the calculations. For example, if 5 folds are used for k-fold CV, and the
#' first 4 are needed to generate these weights, include_folds = c(1,2,3,4)
#' @importFrom stats optim runif
#' @export
#'
#' @examples
#' if (inla_installed()) {
#' spde <- make_mesh(pcod, c("X", "Y"), cutoff = 25)
#'
#' # Set parallel processing if desired:
#' # library(future)
#' # plan(multisession)
#' # depth as quadratic
#' m_cv_1 <- sdmTMB_cv(
#'   density ~ 0 + depth_scaled + depth_scaled2,
#'   data = pcod, spde = spde,
#'   family = tweedie(link = "log"), k_folds = 2
#' )
#' # depth as linear
#' m_cv_2 <- sdmTMB_cv(
#'   density ~ 0 + depth_scaled,
#'   data = pcod, spde = spde,
#'   family = tweedie(link = "log"), k_folds = 2
#' )
#' m_cv_3 <- sdmTMB_cv(
#'   density ~ 1,
#'   data = pcod, spde = spde,
#'   family = tweedie(link = "log"), k_folds = 2
#' )
#' models <- list(m_cv_1, m_cv_2, m_cv_3)
#' weights <- sdmTMB_stacking(model_list = models)
#' }
sdmTMB_stacking <- function(model_list, include_folds = NULL) {
  n_models <- length(model_list)
  if(is.null(include_folds)) {
    n_folds <- max(model_list[[1]]$data$cv_fold)
    include_folds <- seq_len(n_folds)
  }

  # the only quantity we need is the log likelihood
  X <- matrix(0, nrow = nrow(model_list[[1]]$data),
             ncol = n_models)
  for(i in 1:n_models) X[,i] = model_list[[i]]$data$cv_loglik

  # filter out only records in folds of interest
  X <- X[which(model_list[[1]]$data$cv_fold %in% include_folds), ]

  # exponentiate to convert to density in normal space, eqn 5 in Yao et al. (2018)
  X <- exp(X)

  # find model weights that maximize the likelihood
  tot_ll = function(p, X) {
    z <- matrix(exp(p)/sum(exp(p)), ncol=1)
    # see equation below eq 5 in Yao et al. sum over data points of log product of
    #weights and probabilities
    -sum(log(X %*% z)) # neg X because optim will minimize NLL
  }
  # use optimize to find weights that maximize total
  o <- optim(par = runif(n_models), fn = tot_ll, X = X)
  # softmax / normalize to compositions
  weights <- exp(o$par)/sum(exp(o$par))
  return(weights)
}

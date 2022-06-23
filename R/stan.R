#' Extract MCMC samples from a model fit with [tmbstan::tmbstan()].
#'
#' @description
#' `r lifecycle::badge("experimental")`
#'
#' @param object Output from [tmbstan::tmbstan()] run on the `tmb_obj`
#'   element of an [sdmTMB()] model. E.g., `tmbstan(your_model$tmb_obj)`.
#'
#' @return
#' Returns a matrix of parameter samples. Rows correspond to the order
#' of `your_model$tmb_obj$env$last.par.best`. Columns correspond to
#' posterior samples. Is used internally by [predict.sdmTMB()] to make
#' fully Bayesian predictions. See the `tmbstan_model` argument
#' in [predict.sdmTMB()].
#'
#' @examplesIf inla_installed() && require("tmbstan", quietly = TRUE) && ggplot2_installed()
#' \donttest{
#' mesh <- make_mesh(pcod_2011, c("X", "Y"), cutoff = 35) # quite coarse
#'
#' # Fit with marginal maximum likelihood first:
#'
#' fit <- sdmTMB(
#'   density ~ 0 + as.factor(year),
#'   data = pcod_2011, mesh = mesh,
#'   family = tweedie(link = "log"),
#'   priors = sdmTMBpriors(
#'     matern_s = pc_matern(range_gt = 10, sigma_lt = 5),
#'     matern_st = pc_matern(range_gt = 10, sigma_lt = 5),
#'     b = normal(rep(0, 4), scale = rep(10, 4)) # 4 main effects
#'   ),
#'   time = "year"
#' )
#' fit
#'
#' # Create a 'map' vector for TMB to 'fix' kappa at the MLE value
#' # to improve speed of convergence.
#' # Factor NA values cause TMB to fix or map the parameter
#' # at the starting value.
#'
#' pars <- sdmTMB::get_pars(fit)
#' kappa_map <- factor(rep(NA, length(pars$ln_kappa)))
#'
#' # Rebuild model updating some elements:
#' fit_mle <- update(
#'   fit,
#'   control = sdmTMBcontrol(
#'     start = list(
#'       ln_kappa = pars$ln_kappa
#'     ),
#'     map = list(
#'       ln_kappa = kappa_map
#'     )
#'   ),
#'   do_fit = FALSE # no need to actually fit
#' )
#'
#' # Will take a few minutes:
#' library(tmbstan)
#' m_stan <- tmbstan(fit_mle$tmb_obj, iter = 100, chains = 1)
#' print(
#'   m_stan,
#'   pars = c("b_j", "thetaf", "ln_phi", "omega_s[1]", "epsilon_st[1]")
#' )
#'
#' post <- extract_mcmc(m_stan)
#' dim(post)
#'
#' nd <- subset(qcs_grid, year >= 2011)
#' p <- predict(fit_mle, newdata = nd, tmbstan_model = m_stan)
#' p_last <- p[nd$year == max(nd$year), ] # just plot last year
#' pred <- qcs_grid[qcs_grid$year == max(qcs_grid$year), ]
#' pred$est <- apply(exp(p_last), 1, median)
#' pred$lwr <- apply(exp(p_last), 1, quantile, probs = 0.1)
#' pred$upr <- apply(exp(p_last), 1, quantile, probs = 0.9)
#' pred$cv <- apply(exp(p_last), 1, function(x) sd(x) / mean(x))
#'
#' library(ggplot2)
#' ggplot(pred, aes(X, Y, fill = est)) + geom_raster() +
#'   scale_fill_viridis_c(trans = "log")
#' ggplot(pred, aes(X, Y, fill = cv)) + geom_raster() +
#'   scale_fill_viridis_c(trans = "log")
#'
#' index_quantiles <- get_index_sims(p)
#' ggplot(index_quantiles, aes(year, est, ymin = lwr, ymax = upr)) +
#'   geom_line() + geom_ribbon(alpha = 0.5)
#'
#' index_samples <- get_index_sims(p, return_sims = TRUE)
#' ggplot(index_samples, aes(as.factor(year), .value)) +
#'   geom_violin()
#' }
#' @export
extract_mcmc <- function(object) {
  if (!requireNamespace("rstan", quietly = TRUE)) {
    cli_abort("rstan must be installed to use `extract_mcmc()`.")
  }
  post <- rstan::extract(object)
  p_names <- names(post)[-length(names(post))] # exclude "lp__"
  p <- lapply(seq_len(length(post[["lp__"]])), function(i) {
    post_pars <- list()
    for (j in seq_along(p_names)) {
      par_j <- p_names[j]
      if (is.matrix(post[[par_j]])) {
        post_pars[[j]] <- post[[par_j]][i, , drop = TRUE]
      } else {
        post_pars[[j]] <- post[[par_j]][i]
      }
    }
    post_pars
  })
  simplify2array(lapply(p, unlist))
}

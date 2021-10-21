#' Extract MCMC samples from a model fit with [tmbstan::tmbstan()].
#'
#' @description
#' `r lifecycle::badge("experimental")`
#'
#' Returns a matrix of parameter samples. Rows correspond to the order
#' of `your_model$tmb_obj$env$last.par.best`. Columns correspond to
#' posterior samples. Is used internally by [predict.sdmTMB()] to make
#' fully Bayesian predictions. See the `tmbstan_model` argument
#' in [predict.sdmTMB()].
#'
#' @param object Output from [tmbstan::tmbstan()] run on the `tmb_obj`
#'   element of an [sdmTMB()] model. E.g., `tmbstan(your_model$tmb_obj)`.
#' @examples
#'
#' \dontrun{
#' pcod_spde <- make_mesh(pcod, c("X", "Y"), cutoff = 30)
#' plot(pcod_spde)
#'
#' # here we will fix the random field parameters at their approximate
#' # MLEs (maximum likelihood estimates) from a previous fit
#' # to improve speed of convergence:
#' m_tmb <- sdmTMB(density ~ 0 + as.factor(year),
#'   data = pcod, mesh = pcod_spde, family = tweedie(link = "log"), time = "year",
#'   control = sdmTMBcontrol(start = list(ln_kappa = rep(-1.58, 2),
#'     ln_tau_E = -0.15, ln_tau_O = -0.65),
#'     map = list(ln_kappa = rep(factor(NA), 2),
#'       ln_tau_E = factor(NA), ln_tau_O = factor(NA))))
#' m_tmb
#'
#' # will take 3-5 minutes:
#' library(tmbstan)
#' m_stan <- tmbstan(m_tmb$tmb_obj, iter = 200, chains = 1)
#' print(m_stan, pars = c("b_j", "thetaf", "ln_phi", "omega_s[1]", "epsilon_st[1]"))
#'
#' post <- extract_mcmc(m_stan)
#' dim(post)
#'
#' p <- predict(m_tmb, newdata = qcs_grid, tmbstan_model = m_stan)
#' p_last <- p[qcs_grid$year == max(qcs_grid$year), ] # just plot last year
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
    stop("rstan must be installed to use `extract_mcmc()`.", call. = FALSE)
  }
  post <- rstan::extract(object)
  p_names <- names(post)[-length(names(post))] # exclude "lp__"
  p <- lapply(seq_len(nrow(post$b_j)), function(i) {
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

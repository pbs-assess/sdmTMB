#' Extract MCMC samples from a model fit with [tmbstan::tmbstan()].
#'
#' Returns a matrix of parameter samples. Rows correspond to the order
#' of `your_model$tmb_obj$env$last.par.best`. Columns correspond to
#' posterior samples. Can be fed into [predict.sdmTMB()] to make
#' fully Bayesian predictions.
#'
#' @param object Output from [tmbstan::tmbstan()] run on the `tmb_obj`
#'   element of an [sdmTMB()] model. E.g., `tmbstan(your_model$tmb_obj)`.
#' @examples

#' \dontrun{
#' pcod_spde <- make_mesh(pcod, c("X", "Y"), cutoff = 30)
#' plot(pcod_spde)
#'
#' # here we will fix the random field parameters at their approximate
#' # MLEs (maximum likelihood estimates) from a previous fit
#' # to improve speed of convergence:
#' m_tmb <- sdmTMB(density ~ 0 + as.factor(year),
#'   data = pcod, spde = pcod_spde, family = tweedie(link = "log"), time = "year",
#'   start = list(ln_kappa = -1.58, ln_tau_E = -0.15, ln_tau_O = -0.65),
#'   map = list(ln_kappa = factor(NA), ln_tau_E = factor(NA), ln_tau_O = factor(NA)))
#' m_tmb
#'
#' # will take 3-5 minutes:
#' library(tmbstan)
#' m_stan <- tmbstan(m_tmb$tmb_obj, iter = 200, chains = 1)
#' print(m_stan, pars = c("b_j", "thetaf", "ln_phi", "omega_s[1]", "epsilon_st[1]"))
#'
#' post <- extract_mcmc(m_stan)
#' dim(post)
#' }

extract_mcmc <- function(object) {
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

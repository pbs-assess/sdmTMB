qres_tweedie <- function(object, y, mu, ...) {
  p <- stats::plogis(object$model$par[["thetaf"]]) + 1
  dispersion <- exp(object$model$par[["ln_phi"]])

  u <- fishMod::pTweedie(q = y, p = p, mu = mu, phi = dispersion)
  if (p > 1 && p < 2) {
    u[y == 0] <- stats::runif(sum(y == 0), min = 0, max = u[y == 0])
  }
  stats::qnorm(u)
}

qres_binomial <- function(object, y, mu, .n = NULL) {
  # p <- object$family$linkinv(mu) # robust binomial in link space!
  p <- mu
  if (is.null(.n)) .n <- rep(1, length(y))
  mu <- .n * mu
  a <- stats::pbinom(y - 1, .n, p)
  b <- stats::pbinom(y, .n, p)
  u <- stats::runif(n = length(y), min = pmin(a, b), max = pmax(a, b))
  stats::qnorm(u)
}

qres_nbinom2 <- function(object, y, mu, ...) {
  phi <- exp(object$model$par[["ln_phi"]])
  a <- stats::pnbinom(y - 1, size = phi, mu = mu)
  b <- stats::pnbinom(y, size = phi, mu = mu)
  u <- stats::runif(n = length(y), min = a, max = b)
  stats::qnorm(u)
}

# from glmmTMB:
rnbinom1 <- function(n, mu, phi) {
  # var = mu*(1+phi) = mu*(1+mu/k) -> k = mu/phi
  stats::rnbinom(n, mu = mu, size = mu / phi)
}
dnbinom1 <- function(x, mu, phi) {
  stats::dnbinom(x, mu = mu, size = mu / phi)
}
pnbinom1 <- function(q, mu, phi) {
  stats::pnbinom(q, mu = mu, size = mu / phi)
}
qnbinom1 <- function(p, mu, phi) {
  stats::qnbinom(p, mu = mu, size = mu / phi)
}

qres_nbinom1 <- function(object, y, mu, ...) {
  phi <- exp(object$model$par[["ln_phi"]])
  a <- pnbinom1(y - 1, phi = phi, mu = mu)
  b <- pnbinom1(y, phi = phi, mu = mu)
  u <- stats::runif(n = length(y), min = a, max = b)
  stats::qnorm(u)
}

qres_pois <- function(object, y, mu, ...) {
  a <- stats::ppois(y - 1, mu)
  b <- stats::ppois(y, mu)
  u <- stats::runif(n = length(y), min = a, max = b)
  stats::qnorm(u)
}

qres_gamma <- function(object, y, mu, ...) {
  phi <- exp(object$model$par[["ln_phi"]])
  s1 <- phi
  s2 <- mu / s1
  u <- stats::pgamma(q = y, shape = s1, scale = s2)
  stats::qnorm(u)
}

qres_gamma_mix <- function(object, y, mu, ...) {
  p_mix <- plogis(object$model$par[["logit_p_mix"]])
  phi <- exp(object$model$par[["ln_phi"]])
  ratio <- exp(object$model$par[["log_ratio_mix"]])
  s1 <- phi
  s2 <- mu / s1
  s3 <- ratio * s2
  u <- stats::pgamma(q = y, shape = s1, scale = (1-p_mix)*s2 + p_mix*s3)
  stats::qnorm(u)
}

qres_lognormal_mix <- function(object, y, mu, ...) {
  p_mix <- plogis(object$model$par[["logit_p_mix"]])
  dispersion <- exp(object$model$par[["ln_phi"]])
  ratio <- exp(object$model$par[["log_ratio_mix"]])
  u <- stats::plnorm(q = y, meanlog = log((1-p_mix)*mu + p_mix*ratio*mu) - (dispersion^2) / 2, sdlog = dispersion)
  stats::qnorm(u)
}

qres_gaussian <- function(object, y, mu, ...) {
  dispersion <- exp(object$model$par[["ln_phi"]])
  u <- stats::pnorm(q = y, mean = mu, sd = dispersion)
  stats::qnorm(u)
}

qres_lognormal <- function(object, y, mu, ...) {
  dispersion <- exp(object$model$par[["ln_phi"]])
  u <- stats::plnorm(q = y, meanlog = log(mu) - (dispersion^2) / 2, sdlog = dispersion)
  stats::qnorm(u)
}

# https://en.wikipedia.org/wiki/Location%E2%80%93scale_family
pt_ls <- function(q, df, mu, sigma) stats::pt((q - mu) / sigma, df)

qres_student <- function(object, y, mu, ...) {
  dispersion <- exp(object$model$par[["ln_phi"]])
  u <- pt_ls(q = y, df = object$tmb_data$df, mu = mu, sigma = dispersion)
  stats::qnorm(u)
}

qres_beta <- function(object, y, mu, ...) {
  phi <- exp(object$model$par[["ln_phi"]])
  s1 <- mu * phi
  s2 <- (1 - mu) * phi
  u <- stats::pbeta(q = y, shape1 = s1, shape2 = s2)
  stats::qnorm(u)
}

#' Residuals method for sdmTMB models
#'
#' See the residual-checking vignette: `browseVignettes("sdmTMB")` or [on the
#' documentation
#' site](https://pbs-assess.github.io/sdmTMB/articles/residual-checking.html).
#' See notes about types of residuals in 'Details' section below.
#'
#' @param object An [sdmTMB()] model
#' @param type Type of residual. See details.
#' @param mcmc_iter Iterations for MCMC residuals. Will take the last one.
#' @param mcmc_warmup Warmup for MCMC residuals.
#' @param print_stan_model Print the Stan model from MCMC residuals?
#' @param stan_args A list of arguments that will be passed to [rstan::sampling()].
#' @param model Which delta/hurdle model component?
#' @param ... Passed to residual function. Only `n` works for binomial.
#' @export
#' @importFrom stats predict
#' @return A vector of residuals.
#' @seealso [dharma_residuals()]
#' @details
#'
#' Types of residuals currently supported:
#'
#' **`"mle-laplace"`** refers to randomized quantile residuals (Dunn &
#' Smyth 1996), which are also known as probability integral transform (PIT)
#' residuals (Smith 1985). Under model assumptions, these should be distributed
#' as standard normal with the following caveat: the Laplace approximation used
#' for the latent/random effects can cause these residuals to deviate from the
#' standard normal assumption even if the model is consistent with the data
#' (Thygesen et al. 2017). Therefore, **these residuals are fast to calculate
#' but can be unreliable.**
#'
#' **`"mle-mcmc"`** refers to randomized quantile residuals where the fixed
#' effects are fixed at their MLE (maximum likelihoood estimate) values and the
#' random effects are sampled with MCMC via tmbstan/Stan. As proposed in
#' Thygesen et al. (2017) and used in Rufener et al. (2021). Under model
#' assumptions, these should be distributed as standard normal. **These
#' residuals are theoretically preferred over the regular Laplace approximated
#' randomized-quantile residuals, but will be considerably slower to
#' calculate.**
#'
#' Ideally MCMC is run until convergence and then the last iteration can be
#' used for residuals. MCMC samples are defined by `mcmc_iter - mcmc_warmup`.
#' The Stan model can be printed with `print_stan_model = TRUE` to check.
#' The defaults may not be sufficient for many models.
#'
#' **`"mvn-laplace"`** is the same as `"mle-laplace"` except the parameters are
#' based on simulations drawn from the assumed multivariate normal distribution
#' (using the joint precision matrix).
#'
#' **`"response"`** refers to response residuals: observed minus predicted.
#'
#' @references
#' Dunn, P.K. & Smyth, G.K. (1996). Randomized Quantile Residuals. Journal of
#' Computational and Graphical Statistics, 5, 236–244.
#'
#' Smith, J.Q. (1985). Diagnostic checks of non-standard time series models.
#' Journal of Forecasting, 4, 283–291.
#'
#' Rufener, M.-C., Kristensen, K., Nielsen, J.R., and Bastardie, F. 2021.
#' Bridging the gap between commercial fisheries and survey data to model the
#' spatiotemporal dynamics of marine species. Ecological Applications. e02453.
#' \doi{10.1002/eap.2453}
#'
#' Thygesen, U.H., Albertsen, C.M., Berg, C.W., Kristensen, K., and Nielsen, A.
#' 2017. Validation of ecological state space models using the Laplace
#' approximation. Environ Ecol Stat 24(2): 317–339.
#' \doi{10.1007/s10651-017-0372-4}
#'
#' @examples
#' if (inla_installed() &&
#'     require("tmbstan", quietly = TRUE) &&
#'     require("rstan", quietly = TRUE)) {
#'
#'   mesh <- make_mesh(pcod_2011, c("X", "Y"), cutoff = 10)
#'   fit <- sdmTMB(
#'     present ~ as.factor(year) + poly(depth, 3),
#'     data = pcod_2011, mesh = mesh,
#'     family = binomial()
#'   )
#'
#'   # response residuals will be not be normally distributed unless
#'   # the family is Gaussian:
#'   r0 <- residuals(fit, type = "response")
#'   qqnorm(r0)
#'   qqline(r0)
#'
#'   # quick but can have issues because of Laplace approximation:
#'   r1 <- residuals(fit, type = "mle-laplace")
#'   qqnorm(r1)
#'   qqline(r1)
#'
#'   # MCMC-based with fixed effects at MLEs; best but slowest:
#'   set.seed(2938)
#'   r2 <- residuals(fit, type = "mle-mcmc", mcmc_iter = 101, mcmc_warmup = 100)
#'   qqnorm(r2)
#'   qqline(r2)
#'
#'   # Example of passing control arguments to rstan::sampling():
#'   # 11 iterations used for a quick example; don't do this normally
#'   stan_args <- list(control = list(adapt_delta = 0.9, max_treedepth = 12))
#'   r3 <- residuals(
#'     fit, type = "mle-mcmc", mcmc_iter = 11, mcmc_warmup = 10,
#'     stan_args = stan_args
#'   )
#' }

residuals.sdmTMB <- function(object,
                             type = c("mle-laplace", "mle-mcmc", "mvn-laplace", "response", "pearson"),
                             mcmc_iter = 500, mcmc_warmup = 250,
                             print_stan_model = FALSE,
                             stan_args = NULL,
                             model = c(1, 2),
                             ...) {

  model <- as.integer(model[[1]])
  if ("visreg_model" %in% names(object)) {
    model <- object$visreg_model
  }

  # retrieve function that called this:
  sys_calls <- unlist(lapply(sys.calls(), deparse))
  visreg_call <- any(grepl("setupV", substr(sys_calls, 1, 7)))

  type <- match.arg(type)
  if (!visreg_call) {
    msg <- c(
      "We recommend using the slower `type = 'mle-mcmc'` for final inference.",
      "See the ?residuals.sdmTMB 'Details' section."
    )
    if (type == "randomized-quantile") cli_inform(msg)
  }

  fam <- object$family$family
  nd <- NULL
  est_column <- "est"
  linkinv <- object$family$linkinv
  if (isTRUE(object$family$delta)) {
    fam <- fam[[model]]
    linkinv <-  object$family[[model]]$linkinv
    nd <- object$data
    est_column <- if (model == 1L) "est1" else "est2"
  }
  res_func <- switch(fam,
    gaussian = qres_gaussian,
    binomial = qres_binomial,
    tweedie  = qres_tweedie,
    Beta     = qres_beta,
    Gamma    = qres_gamma,
    nbinom2  = qres_nbinom2,
    nbinom1  = qres_nbinom1,
    poisson  = qres_pois,
    student  = qres_student,
    lognormal  = qres_lognormal,
    gamma_mix = qres_gamma_mix,
    lognormal_mix = qres_lognormal_mix,
    cli_abort(paste(fam, "not yet supported."))
  )

  if (type %in% c("mle-laplace", "response", "pearson")) {
    mu <- tryCatch({linkinv(predict(object, newdata = nd)[[est_column]])}, # newdata = NULL; fast
      error = function(e) NA)
    if (is.na(mu[[1]])) {
      mu <- linkinv(predict(object, newdata = object$data)[[est_column]]) # not newdata = NULL
    }
  } else if (type == "mvn-laplace") {
    mu <- linkinv(predict(object, nsim = 1L, model = model)[, 1L, drop = TRUE])
  } else if (type == "mle-mcmc") {
    if (!requireNamespace("rstan", quietly = TRUE))
      cli_abort("rstan must be installed.")
    if (!requireNamespace("tmbstan", quietly = TRUE))
      cli_abort("tmbstan must be installed.")

    assert_that(is.numeric(mcmc_iter))
    assert_that(is.numeric(mcmc_warmup))
    assert_that(mcmc_warmup < mcmc_iter)
    # from https://github.com/mcruf/LGNB/blob/8aba1ee2df045c2eb45e124d5a753e8f1c6e865a/R/Validation_and_Residuals.R
    # get names of random effects in the model
    obj <- object$tmb_obj
    random <- unique(names(obj$env$par[obj$env$random]))
    if (isTRUE(object$reml)) {
      msg <- c(
        "Please refit your model with `reml = FALSE` to use MCMC-MLE residuals.",
        "You can use `fit_ml <- update(your_reml_fit, reml = FALSE)` to do this."
      )
      cli_abort(msg)
    }
    # get (logical) non random effects indices:
    pl <- as.list(object$sd_report, "Estimate")
    fixed <- !(names(pl) %in% random)
    # fix non-random parameters to their estimated values:
    map <- lapply(pl[fixed], function(x) factor(rep(NA, length(x))))
    # construct corresponding new function object:
    obj <- TMB::MakeADFun(obj$env$data, pl, map = map, DLL = "sdmTMB")
    # run MCMC to get posterior sample of random effects given data:

    args <- list(obj = obj, chains = 1L, iter = mcmc_iter, warmup = mcmc_warmup)
    args <- c(args, stan_args)

    samp <- do.call(tmbstan::tmbstan, args)
    if (print_stan_model) print(samp)
    obj_mle <- object
    obj_mle$tmb_obj <- obj
    obj_mle$tmb_map <- map
    pred <- predict(obj_mle, tmbstan_model = samp, model = model, nsim = 1L) # only use last
    # pred <- as.numeric(pred[, ncol(pred), drop = TRUE])
    mu <- linkinv(pred)
  } else {
    cli_abort("residual type not implemented")
  }

  y <- object$response
  y <- y[, model, drop = TRUE] # in case delta
  # e.g., visreg, prediction has already removed NA mu:
  if (sum(is.na(y)) > 0 && length(mu) < length(y)) {
    y <- y[!is.na(y)]
  }

  # for binomial proportion with weights = N:
  size <- object$tmb_data$size
  prop_binomial <- !all(size == 1)

  if (type == "response") {
    if (!prop_binomial) r <- y - mu else r <- y / size - mu
  } else if (type == "mle-laplace" || type == "mvn-laplace") {
    r <- res_func(object, y, mu, .n = size, ...)
  } else if (type == "mle-mcmc") {
    r <- res_func(obj_mle, y, mu, .n = size, ...)
  } else if (type == "pearson") {
    if (is.null(v <- family(object)$variance)) {
      cli_abort(c("Variance function undefined for family;",
        "cannot compute Pearson residuals"))
    }
    # FIXME: add sigma function, for now just binomial
    if (length(formals(v)) > 1)
      cli_abort("sdmTMB currently only supports variance functions with 1 argument")
    # vv <- switch(length(formals(v)),
    #   v(fitted(object)),
    #   v(fitted(object), sigma(object)),
    #   stop("variance function should take 1 or 2 arguments"))
    vv <- v(mu)
    wts <- if (prop_binomial) size else object$tmb_data$weights_i
    if (!prop_binomial) r <- y - mu else r <- y / size - mu
    r <- r / sqrt(vv)
    if (!is.null(wts)) r <- r * sqrt(wts)
  } else {
    cli_abort("residual type not implemented")
  }
  r
}

# from:
# https://bbolker.github.io/mixedmodels-misc/glmmFAQ.html#testing-for-overdispersioncomputing-overdispersion-factor
check_overdisp <- function(object) {
  rdf <- stats::df.residual(object)
  rp <- stats::residuals(object, type = "pearson")
  pearson_chisq <- sum(rp^2)
  prat <- pearson_chisq / rdf
  pval <- stats::pchisq(pearson_chisq, df = rdf, lower.tail = FALSE)
  data.frame(chisq = pearson_chisq, ratio = prat, rdf = rdf, p = pval)
}

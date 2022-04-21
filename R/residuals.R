qres_tweedie <- function(object, y, mu) {
  p <- stats::plogis(object$model$par[["thetaf"]]) + 1
  dispersion <- exp(object$model$par[["ln_phi"]])

  u <- fishMod::pTweedie(q = y, p = p, mu = mu, phi = dispersion)
  if (p > 1 && p < 2)
    u[y == 0] <- stats::runif(sum(y == 0), min = 0, max = u[y == 0])
  stats::qnorm(u)
}

qres_binomial <- function(object, y, mu, n = NULL) {
  p <- object$family$linkinv(mu) # robust binomial in link space!
  if (is.null(n)) n <- rep(1, length(y))
  y <- n * y
  a <- stats::pbinom(y - 1, n, p)
  b <- stats::pbinom(y, n, p)
  u <- stats::runif(n = length(y), min = pmin(a, b), max = pmax(a, b))
  stats::qnorm(u)
}

qres_nbinom2 <- function(object, y, mu) {
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

qres_nbinom1 <- function(object, y, mu) {
  phi <- exp(object$model$par[["ln_phi"]])
  a <- pnbinom1(y - 1, phi = phi, mu = mu)
  b <- pnbinom1(y, phi = phi, mu = mu)
  u <- stats::runif(n = length(y), min = a, max = b)
  stats::qnorm(u)
}

qres_pois <- function(object, y, mu) {
  a <- stats::ppois(y - 1, mu)
  b <- stats::ppois(y, mu)
  u <- stats::runif(n = length(y), min = a, max = b)
  stats::qnorm(u)
}

qres_gamma <- function(object, y, mu) {
  phi <- exp(object$model$par[["ln_phi"]])
  s1 <- phi
  s2 <- mu / s1
  u <- stats::pgamma(q = y, shape = s1, scale = s2)
  stats::qnorm(u)
}

qres_gaussian <- function(object, y, mu) {
  dispersion <- exp(object$model$par[["ln_phi"]])
  u <- stats::pnorm(q = y, mean = mu, sd = dispersion)
  stats::qnorm(u)
}

qres_lognormal <- function(object, y, mu) {
  dispersion <- exp(object$model$par[["ln_phi"]])
  u <- stats::plnorm(q = y, meanlog = log(mu) - (dispersion^2)/2, sdlog = dispersion)
  stats::qnorm(u)
}

# https://en.wikipedia.org/wiki/Location%E2%80%93scale_family
pt_ls <- function(q, df, mu, sigma) stats::pt((q - mu)/sigma, df)

qres_student <- function(object, y, mu) {
  dispersion <- exp(object$model$par[["ln_phi"]])
  u <- pt_ls(q = y, df = object$tmb_data$df, mu = mu, sigma = dispersion)
  stats::qnorm(u)
}

qres_beta <- function(object, y, mu) {
  phi <- exp(object$model$par[["ln_phi"]])
  s1 <- mu * phi
  s2 <- (1 - mu) * phi
  u <- stats::pbeta(q = y, shape1 = s1, shape2 = s2)
  stats::qnorm(u)
}

#' Residuals method for sdmTMB models
#'
#' These residuals are randomized quantile residuals (Dunn & Smyth 1996), which
#' are also known as probability integral transform (PIT) residuals (Smith
#' 1985). See the residual-checking vignette: `browseVignettes("sdmTMB")` or [on
#' the documentation
#' site](https://pbs-assess.github.io/sdmTMB/articles/residual-checking.html).
#'
#' @param object An [sdmTMB()] model
#' @param type Type of residual. Residual at the MLE or based on simulations
#'   from the joint precision matrix are available.
#' @param ... Passed to residual function. Only `n` works for binomial.
#' @export
#' @importFrom stats predict
#' @return A vector of randomized quantile residuals.
#' @seealso [dharma_residuals()]
#' @references
#' Dunn, P.K. & Smyth, G.K. (1996). Randomized Quantile Residuals. Journal of
#' Computational and Graphical Statistics, 5, 236–244.
#'
#' Smith, J.Q. (1985). Diagnostic checks of non-standard time series models.
#' Journal of Forecasting, 4, 283–291.
residuals.sdmTMB <- function(object, type = c("mle", "sim"), ...) {
  if (isTRUE(object$family$delta)) {
    nice_stop("`residuals.sdmTMB()` is not setup to work with delta models yet. ",
      "Try `dharma_residuals()`.")
  }
  # message("Consider using `dharma_residuals()` instead.")
  type <- match.arg(type)
  res_func <- switch(object$family$family,
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
    nice_stop(paste(object$family$family, "not yet supported."))
  )
  if (type == "mle") {
    mu <- object$family$linkinv(predict(object, newdata = NULL)$est)
  } else if (type == "sim") {
    mu <- object$family$linkinv(predict(object, nsim = 1L)[, 1L, drop = TRUE])
  } else {
    nice_stop("`type` not implemented")
  }
  y <- object$response
  y <- y[,1,drop=TRUE]
  res_func(object, y, mu, ...)
}

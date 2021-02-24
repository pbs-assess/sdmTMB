qres_tweedie <- function(object, y, mu) {
  p <- stats::plogis(object$model$par[["thetaf"]]) + 1
  dispersion <- exp(object$model$par[["ln_phi"]])
  u <- fishMod::pTweedie(q = y, p = p, mu = mu, phi = dispersion)
  if (p > 1 && p < 2)
    u[y == 0] <- stats::runif(sum(y == 0), min = 0, max = u[y == 0])
  stats::qnorm(u)
}

qres_binomial <- function(object, y, mu, n=NULL) {
  if(is.null(n)) n <- rep(1, length(y))
  a <- stats::pbinom(n-y, n, mu)
  b <- stats::pbinom(y, n, mu)
  u <- stats::runif(n = length(y), min = min(a,b), max = max(a,b))
  stats::qnorm(u)
}

qres_nbinom2 <- function(object, y, mu) {
  phi <- exp(object$model$par[["ln_phi"]])
  a <- stats::pnbinom(y - 1, size = phi, mu = mu)
  b <- stats::pnbinom(y, size = phi, mu = mu)
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

#' @export
#' @importFrom stats predict
residuals.sdmTMB <- function(object, ...) {
  res_func <- switch(object$family$family,
    gaussian = qres_gaussian,
    binomial = qres_binomial,
    tweedie  = qres_tweedie,
    Beta     = qres_beta,
    Gamma    = qres_gamma,
    nbinom2  = qres_nbinom2,
    poisson  = qres_pois,
    student  = qres_student
  )
  y <- object$response
  mu <- object$family$linkinv(predict(object)$est)
  res_func(object, y, mu, ...)
}

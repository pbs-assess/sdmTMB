qres_tweedie <- function (object, y, mu) {
  w <- 1 # weights
  p <- stats::plogis(object$model$par[["thetaf"]]) + 1
  dispersion <- exp(object$model$par[["ln_phi"]])
  u <- tweedie::ptweedie(q = y, power = p, mu = mu, phi = dispersion/w)
  if (p > 1 && p < 2)
    u[y == 0] <- stats::runif(sum(y == 0), min = 0, max = u[y == 0])
  stats::qnorm(u)
}

qres_binomial <- function (object, y, mu) {
  n <- rep(1, length(y))
  a <- stats::pbinom(y - 1, n, mu)
  b <- stats::pbinom(y, n, mu)
  u <- stats::runif(n = length(y), min = a, max = b)
  stats::qnorm(u)
}

qres_gamma <- function (object, y, mu) {
  stop("Not finished.")
  # df <- stats::df.residual(object)
  # w <- 1
  # # dispersion <- sum(w * ((y - mu)/mu)^2)/df
  # dispersion <- 1/exp(object$fit$par[["betad"]])
  # logp <- stats::pgamma((w * y)/mu/dispersion, w/dispersion, log.p = TRUE)
  # stats::qnorm(logp, log.p = TRUE)
}

qres_gaussian <- function (object, y, mu) {
  dispersion <- exp(object$model$par[["ln_phi"]])
  u <- stats::pnorm(q = y, mean = mu, sd = dispersion)
  stats::qnorm(u)
}

#' @export
residuals.sdmTMB <- function(object, ...) {
  res_func <- switch(object$family$family,
    gaussian = qres_gaussian,
    binomial = qres_binomial,
    tweedie  = qres_tweedie
  )
  y <- object$response
  mu <- object$family$linkinv(predict(object)$est)
  res_func(object, y, mu)
}

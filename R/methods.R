#' @import methods
#' @export
summary.sdmTMB <- function(object, ..., digits) {
  print(object, ...)
}

mround <- function(x, digits) {
  sprintf(paste0("%.", digits, "f"), round(x, digits))
}

#' Extract the number of observations of an sdmTMB model
#'
#' @param object The fitted sdmTMB model object
#' @importFrom stats nobs
#' @export
#' @noRd
nobs.sdmTMB <- function(object, ...) {
    sum(!is.na(object$data[all.vars(object$formula[[1]])[1]]))
}

#' Get fitted values from an sdmTMB model
#'
#' @param object The fitted sdmTMB model object
#' @importFrom stats predict
#' @export
#' @noRd
fitted.sdmTMB <- function(object, ...) {
  predict(object, type = "response")$est
}

#' Extract the log likelihood of a sdmTMB model
#'
#' @param object The fitted sdmTMB model object
#' @importFrom stats logLik
#' @export
#' @noRd
logLik.sdmTMB <- function(object, ...) {
  val <- -object$model$objective
  nobs <- nobs.sdmTMB(object)
  df <- length(object$model$par) # fixed effects only
  if (isTRUE(object$reml)) {
    s <- as.list(object$sd_report, "Estimate")
    df <- df + length(s$b_j) + length(s$b_j2)
  }
  structure(val,
    nobs = nobs, nall = nobs, df = df,
    class = "logLik"
  )
}

#' Extract the AIC of a sdmTMB model
#'
#' @param fit The fitted sdmTMB model
#' @param scale The scale (note used)
#' @param k Penalization parameter, defaults to 2
#' @param ... Anything else
#' @noRd
#'
#' @export
extractAIC.sdmTMB <- function(fit, scale, k = 2, ...) {
  L <- logLik(fit)
  edf <- attr(L, "df")
  return(c(edf, c(-2 * L + k * edf)))
}

#' @importFrom stats family
#' @export
family.sdmTMB <- function (object, ...) {
  if ("visreg_model" %in% names(object)) {
    return(object$family[[object$visreg_model]])
  } else {
    return(object$family)
  }
}

##' @importFrom nlme fixef
#' @method fixef sdmTMB
#' @export
fixef.sdmTMB <- function(object, ...) {
  .t <- tidy(object)
  bhat <- .t$estimate
  names(bhat) <- .t$term
  bhat
}

##' @importFrom nlme ranef
#' @method ranef sdmTMB
#' @export
ranef.sdmTMB <- function(object, ...) {
  .t <- tidy(object, "ran_vals", conf.int = FALSE)
  terms <- unlist(lapply(strsplit(.t$term,"_"), getElement, 1))
  est <- .t$estimate
  cond <- list()
  for(i in seq_along(unique(terms))) {
    cond[[unique(terms)[i]]] =
      data.frame("Intercept" = est[which(terms == unique(terms)[i])], stringsAsFactors = FALSE)
  }
  return(list(cond = cond))
}

#' @importFrom stats df.residual
#' @method df.residual sdmTMB
#' @export
df.residual.sdmTMB <- function(object, ...) {
  nobs(object) - length(object$model$par)
}

#' @export
formula.sdmTMB <- function (x, ...) {
  if (length(x$formula) > 1L) {
    if ("visreg_model" %in% names(x)) {
      return(x$formula[[x$visreg_model]])
    } else {
      return(x$formula)
    }
  } else {
    return(x$formula[[1]])
  }
}

#' @importFrom stats vcov
#' @export
vcov.sdmTMB <- function(object, ...) {
  vc <- object$sd_report$cov.fixed
  rn <- rownames(vc)
  bj <- grepl("^b_j", rn)
  vc <- vc[bj, bj]
  b <- tidy(object)
  stopifnot(nrow(b) == nrow(vc))
  rownames(vc) <- b$term
  colnames(vc) <- b$term
  vc
}

#' @importFrom stats terms
#' @export
terms.sdmTMB <- function(x, ...) {
  # DELTA FIXME: hardcoded to model 1!
  class(x) <- "glm" # fake
  out <- stats::terms(x)
  out[[1]]
}

#' Calculate effects
#'
#' Used by effects package
#'
#' @inheritParams effects::Effect
#'
#' @importFrom stats formula poisson
#'
#' @return
#' Output from [effects::effect()]. Can then be plotted with with associated
#' `plot()` method.
#'
#' @rawNamespace if(getRversion() >= "3.6.0") {
#'   S3method(effects::Effect, sdmTMB)
#' } else {
#'   export(Effect.sdmTMB)
#' }
#' @examplesIf inla_installed() && require("effects", quietly = TRUE)
#' fit <- sdmTMB(present ~ depth_scaled, data = pcod_2011, family = binomial(),
#'   spatial = "off")
#' effects::effect("depth_scaled", fit)
#' plot(effects::effect("depth_scaled", fit))
Effect.sdmTMB <- function(focal.predictors, mod, ...) {
  if (!requireNamespace("effects", quietly = TRUE)) {
    cli_abort("Please install the effects package")
  }

  vc <- vcov(mod)
  b <- tidy(mod)

  dummyfuns <- list(
    variance = function(mu) mu,
    initialize = expression(mustart = y + 0.1),
    dev.resids = function(...) stats::poisson()$dev.res(...)
  )

  fam <- family(mod)

  # from glmmTMB:
  dummyfuns <- list(
    variance = function(mu) mu,
    initialize = expression(mustart <- y + 0.1),
    dev.resids = function(...) poisson()$dev.res(...)
  )
  for (i in names(dummyfuns)) {
    if (is.null(fam[[i]])) fam[[i]] <- dummyfuns[[i]]
  }

  coefs <- b$estimate
  names(coefs) <- b$term
  args <- list(
    call = mod$call,
    coefficients = coefs,
    vcov = vc,
    family = fam,
    formula = formula(mod)
  )
  effects::Effect.default(focal.predictors, mod, ..., sources = args)
}

# get_term_names <- function(model) {
#   .names <- gsub(" ", "", labels(terms(model)))
#   .names
# }

#' @export
model.frame.sdmTMB <- function(formula, ...) {
  as.data.frame(formula$data) # no tibbles!
}

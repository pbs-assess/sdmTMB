#' Estimated marginal means with the \pkg{emmeans} package with \pkg{sdmTMB}
#'
#' @description
#' Methods for using the \pkg{emmeans} package with \pkg{sdmTMB}. The
#' \pkg{emmeans} package computes estimated marginal means for the fixed
#' effects.
#'
#' @name emmeans.sdmTMB
#'
#' @references
#' \url{https://aosmith.rbind.io/2019/03/25/getting-started-with-emmeans/}
#'
#' @examplesIf inla_installed() && require("emmeans", quietly = TRUE)
#' mesh <- make_mesh(pcod_2011, c("X", "Y"), cutoff = 20)
#' fit <- sdmTMB(
#'   present ~ as.factor(year),
#'   data = pcod_2011, mesh = mesh,
#'   family = binomial()
#' )
#' fit
#' emmeans::emmeans(fit, ~ year)
#' emmeans::emmeans(fit, pairwise ~ year)
#' emmeans::emmeans(fit, pairwise ~ year, type = "response")
#' emmeans::emmeans(fit, pairwise ~ year, adjust = "none")
#'
#' e <- emmeans::emmeans(fit, ~ year)
#' plot(e)
#'
#' e <- emmeans::emmeans(fit, pairwise ~ year)
#' confint(e)
#' summary(e, infer = TRUE)
#' as.data.frame(e)
#'
#' # interaction of factor with continuous predictor:
#' fit2 <- sdmTMB(
#'   present ~ depth_scaled * as.factor(year),
#'   data = pcod_2011, mesh = mesh,
#'   family = binomial()
#' )
#' fit2
#' # slopes for each level:
#' emmeans::emtrends(fit2, ~ year, var = "depth_scaled")
#' # test difference in slopes:
#' emmeans::emtrends(fit2, pairwise ~ year, var = "depth_scaled")
#' emmeans::emmip(fit2, year ~ depth_scaled,
#'   at = list(depth_scaled = seq(-2.5, 2.5, length.out = 50)), CIs = TRUE)

NULL # don't document functions below

recover_data.sdmTMB <- function(object, ...) {
  fcall <- stats::getCall(object)
  if (!requireNamespace("emmeans", quietly = TRUE)) {
    cli_abort("Please install the emmeans package to use this function")
  }
  emmeans::recover_data(
    fcall,
    stats::delete.response(terms(object)),
    attr(model.frame(object), "na.action"), ...
  )
}

# with help from emm_basis.glmmTMB
emm_basis.sdmTMB <- function(object, trms, xlev, grid, ...) {
  V <- vcov(object)
  misc <- list()
  fam <- family(object)
  misc <- emmeans::.std.link.labels(fam, misc)
  contrasts <- attr(model.matrix(object), "contrasts")
  contrasts <- contrasts[names(contrasts) %in% all.vars(terms(object))]
  m <- model.frame(trms, grid, na.action = stats::na.pass, xlev = xlev)
  X <- model.matrix(trms, m, contrasts.arg = contrasts)
  bhat <- fixef(object)
  if (length(bhat) < ncol(X)) {
    kept <- match(names(bhat), dimnames(X)[[2]])
    bhat <- NA * X[1, ]
    bhat[kept] <- fixef(object)
    modmat <- model.matrix(
      trms, model.frame(object),
      contrasts.arg = contrasts
    )
    nbasis <- estimability::nonest.basis(modmat)
  } else {
    nbasis <- estimability::all.estble
  }
  dfargs <- list(df = df.residual(object))
  dffun <- function(k, dfargs) dfargs$df
  named_list(X, bhat, nbasis, V, dffun, dfargs, misc)
}

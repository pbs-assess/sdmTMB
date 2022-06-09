# modified from glmmTMB
# extra stuff for Effects package, class, etc.
add_to_family <- function(x) {
  # x <- c(x, list(link = link), make.link(link))
  # Effect.default/glm.fit
  if (is.null(x$aic)) {
    x <- c(x, list(aic = function(...) NA_real_))
  }
  if (is.null(x$initialize)) {
    x <- c(x, list(initialize = expression({
      mustart <- y + 0.1
    })))
  }
  if (is.null(x$dev.resids)) {
    # can't return NA, glm.fit is unhappy
    x <- c(x, list(dev.resids = function(y, mu, wt) {
      rep(0, length(y))
    }))
  }
  class(x) <- "family"
  x
}

#' Additional families
#'
#' Additional families compatible with [sdmTMB()].
#'
#' @param link Link.
#' @export
#' @rdname families
#' @name Families
#'
#' @return
#' A list with elements `family`, `link`, `linkfun`, and `linkinv`.
#'
#' @examples
#' Beta(link = "logit")
Beta <- function(link = "logit") {
  linktemp <- substitute(link)
  if (!is.character(linktemp))
    linktemp <- deparse(linktemp)
  okLinks <- c("logit")
  if (linktemp %in% okLinks)
    stats <- stats::make.link(linktemp)
  else if (is.character(link))
    stats <- stats::make.link(link)
  x <- c(list(family = "Beta", link = linktemp), stats)
  add_to_family(x)
}

#' @export
#' @rdname families
#' @examples
#' lognormal(link = "log")
lognormal <- function(link = "log") {
  linktemp <- substitute(link)
  if (!is.character(linktemp))
    linktemp <- deparse(linktemp)
  okLinks <- c("identity", "log", "inverse")
  if (linktemp %in% okLinks)
    stats <- stats::make.link(linktemp)
  else if (is.character(link))
    stats <- stats::make.link(link)
  x <- c(list(family = "lognormal", link = linktemp), stats)
  add_to_family(x)
}

#' @details
#' The `nbinom2` negative binomial parameterization is the NB2 where the variance grows
#' quadratically with the mean (Hilbe 2011).
#' @references
#' Hilbe, J. M. (2011). Negative binomial regression. Cambridge University Press.
#' @export
#' @examples
#' nbinom2(link = "log")
#' @rdname families
nbinom2 <- function(link = "log") {
  linktemp <- substitute(link)
  if (!is.character(linktemp))
    linktemp <- deparse(linktemp)
  okLinks <- c("log")
  if (linktemp %in% okLinks)
    stats <- stats::make.link(linktemp)
  else if (is.character(link))
    stats <- stats::make.link(link)
  x <- c(list(family = "nbinom2", link = linktemp), stats)
  add_to_family(x)
}

#' @details
#' The `nbinom1` negative binomial parameterization lets the variance grow
#' linearly with the mean (Hilbe 2011).
#' @export
#' @examples
#' nbinom1(link = "log")
#' @rdname families
nbinom1 <- function(link = "log") {
  linktemp <- substitute(link)
  if (!is.character(linktemp))
    linktemp <- deparse(linktemp)
  okLinks <- c("log")
  if (linktemp %in% okLinks)
    stats <- stats::make.link(linktemp)
  else if (is.character(link))
    stats <- stats::make.link(link)
  x <- c(list(family = "nbinom1", link = linktemp), stats)
  add_to_family(x)
}

#' @export
#' @examples
#' truncated_nbinom2(link = "log")
#' @rdname families
truncated_nbinom2 <- function(link = "log") {
  linktemp <- substitute(link)
  if (!is.character(linktemp))
    linktemp <- deparse(linktemp)
  okLinks <- c("log")
  if (linktemp %in% okLinks)
    stats <- stats::make.link(linktemp)
  else if (is.character(link))
    stats <- stats::make.link(link)

  list(family = "truncated_nbinom2", link = linktemp, linkfun = stats$linkfun,
    linkinv = stats$linkinv)
}

#' @export
#' @examples
#' truncated_nbinom1(link = "log")
#' @rdname families
truncated_nbinom1 <- function(link = "log") {
  linktemp <- substitute(link)
  if (!is.character(linktemp))
    linktemp <- deparse(linktemp)
  okLinks <- c("log")
  if (linktemp %in% okLinks)
    stats <- stats::make.link(linktemp)
  else if (is.character(link))
    stats <- stats::make.link(link)

  list(family = "truncated_nbinom1", link = linktemp, linkfun = stats$linkfun,
    linkinv = stats$linkinv)
}

#' @param df Student-t degrees of freedom fixed value parameter.
#' @export
#' @details
#' The degrees of freedom parameter is currently not estimated and is fixed at `df`.
#' @rdname families
#' @examples
#' student(link = "identity")
student <- function(link = "identity", df = 3) {
  linktemp <- substitute(link)
  if (!is.character(linktemp))
    linktemp <- deparse(linktemp)
  okLinks <- c("identity", "log", "inverse")
  if (linktemp %in% okLinks)
    stats <- stats::make.link(linktemp)
  else if (is.character(link))
    stats <- stats::make.link(link)

  x <- c(list(family = "student", link = linktemp, df = df), stats)
  add_to_family(x)
}

#' @export
#' @examples
#' tweedie(link = "log")
#' @rdname families
tweedie <- function(link = "log") {
  linktemp <- substitute(link)
  if (!is.character(linktemp))
    linktemp <- deparse(linktemp)
  okLinks <- c("log", "identity")
  if (linktemp %in% okLinks)
    stats <- stats::make.link(linktemp)
  else if (is.character(link))
    stats <- stats::make.link(link)

  x <- c(list(family = "tweedie", link = linktemp), stats)
  add_to_family(x)
}

#' @export
#' @examples
#' censored_poisson(link = "log")
#' @rdname families
censored_poisson <- function(link = "log") {
  linktemp <- substitute(link)
  if (!is.character(linktemp))
    linktemp <- deparse(linktemp)
  okLinks <- c("log")
  if (linktemp %in% okLinks)
    stats <- stats::make.link(linktemp)
  else if (is.character(link))
    stats <- stats::make.link(link)

  list(family = "censored_poisson", link = linktemp, linkfun = stats$linkfun,
    linkinv = stats$linkinv)
}

#' @param link1 Link for first part of delta/hurdle model.
#' @param link2 Link for second part of delta/hurdle model.
#' @export
#' @importFrom stats Gamma binomial
#' @examples
#' delta_gamma()
#' @rdname families
delta_gamma <- function(link1 = "logit", link2 = "log") {
  link1 <- match.arg(link1)
  link2 <- match.arg(link2)
  f1 <- binomial(link = "logit")
  f2 <- Gamma(link = "log")
  list(f1, f2, delta = TRUE, link = c("logit", "log"),
    family = c("binomial", "Gamma"),
    clean_name = "delta_gamma(link1 = 'logit', link2 = 'log')")
}

#' @export
#' @examples
#' delta_lognormal()
#' @rdname families
delta_lognormal <- function(link1 = "logit", link2 = "log") {
  link1 <- match.arg(link1)
  link2 <- match.arg(link2)
  f1 <- binomial(link = "logit")
  f2 <- lognormal(link = "log")
  list(f1, f2, delta = TRUE, link = c("logit", "log"),
    family = c("binomial", "lognormal"),
    clean_name = "delta_lognormal(link1 = 'logit', link2 = 'log')")
}

#' @export
#' @examples
#' delta_truncated_nbinom2()
#' @rdname families
delta_truncated_nbinom2 <- function(link1 = "logit", link2 = "log") {
  link1 <- match.arg(link1)
  link2 <- match.arg(link2)
  f1 <- binomial(link = "logit")
  f2 <- truncated_nbinom2(link = "log")
  list(f1, f2, delta = TRUE, link = c("logit", "log"),
    family = c("binomial", "truncated_nbinom2"),
    clean_name = "delta_truncated_nbinom2(link1 = 'logit', link2 = 'log')")
}

#' @export
#' @examples
#' delta_truncated_nbinom1()
#' @rdname families
delta_truncated_nbinom1 <- function(link1 = "logit", link2 = "log") {
  link1 <- match.arg(link1)
  link2 <- match.arg(link2)
  f1 <- binomial(link = "logit")
  f2 <- truncated_nbinom1(link = "log")
  list(f1, f2, delta = TRUE, link = c("logit", "log"),
    family = c("binomial", "truncated_nbinom1"),
    clean_name = "delta_truncated_nbinom1(link1 = 'logit', link2 = 'log')")
}

# @examples
# delta_poisson_link_gamma()
# @rdname families
# @details `delta_poisson_link_gamma()` is the Poisson-link (complementary
#   log-log) delta model (Thorson 2018).
# @references
# Thorson, J. T. (2018). Three problems with the conventional delta-model for
# biomass sampling data, and a computationally efficient alternative. Canadian
# Journal of Fisheries and Aquatic Sciences, 75(9), 1369-1382.
# \doi{10.1139/cjfas-2017-0266}
# @export
# delta_poisson_link_gamma <- function(link1 = "log", link2 = "log") {
#   cli_abort("`delta_poisson_link_gamma()` is experimental and not all functions work with it")
#   cli_inform("Index calculations may not be correct with the `delta_poisson_link_gamma()` family yet")
#   link1 <- match.arg(link1)
#   link2 <- match.arg(link2)
#   f1 <- binomial(link = "log")
#   f2 <- truncated_nbinom1(link = "log")
#   list(f1, f2, delta = TRUE, link = c("log", "log"),
#     family = c("binomial", "Gamma"), type = "poisson_link_delta",
#     clean_name = "delta_poisson_link_gamma(link1 = 'logit', link2 = 'log')")
# }

#' Additional families
#'
#' Additional families compatible with [sdmTMB()].
#'
#' @param link The link.
#' @export
#' @rdname families
#' @name Families
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
  list(family = "Beta", link = linktemp, linkfun = stats$linkfun,
    linkinv = stats$linkinv)
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

  list(family = "lognormal", link = linktemp, linkfun = stats$linkfun,
    linkinv = stats$linkinv)
}

#' @details
#' The negative binomial parameterization is the NB2 where the variance grows
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

  list(family = "nbinom2", link = linktemp, linkfun = stats$linkfun,
    linkinv = stats$linkinv)
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

  list(family = "student", link = linktemp, linkfun = stats$linkfun,
    linkinv = stats$linkinv, df = df)
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

  list(family = "tweedie", link = linktemp, linkfun = stats$linkfun,
    linkinv = stats$linkinv)
}

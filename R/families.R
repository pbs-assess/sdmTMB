#' Additional families
#'
#' Additional families compatible with [sdmTMB()].
#'
#' @param link The link.
#' @export
#' @examples
#' tweedie(link = "log")
#' @rdname families
tweedie <- function(link = "log") {
  linktemp <- substitute(link)
  if (!is.character(linktemp))
    linktemp <- deparse(linktemp)
  okLinks <- c("inverse", "log", "identity")
  if (linktemp %in% okLinks)
    stats <- stats::make.link(linktemp)
  else if (is.character(link))
    stats <- stats::make.link(link)

  list(family = "tweedie", link = linktemp, linkfun = stats$linkfun,
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

#' @export
#' @details
#' The degrees of freedom parameter for the Student-t distribution is currently
#' fixed at 3.
#' @rdname families
#' @examples
#' student(link = "identity")
student <- function(link = "identity") {
  linktemp <- substitute(link)
  if (!is.character(linktemp))
    linktemp <- deparse(linktemp)
  okLinks <- c("identity", "log", "inverse")
  if (linktemp %in% okLinks)
    stats <- stats::make.link(linktemp)
  else if (is.character(link))
    stats <- stats::make.link(link)

  list(family = "student", link = linktemp, linkfun = stats$linkfun,
    linkinv = stats$linkinv)
}

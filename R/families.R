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
  list(family = "nbinom2", link = linktemp, linkfun = stats$linkfun,
    linkinv = stats$linkinv)
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
  list(family = "nbinom1", link = linktemp, linkfun = stats$linkfun,
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

.make_delta_family <- function(link1 = "logit", link2 = "log",
  ok_links = c("logit", "log")) {

}

.make_delta <- function(link1, link2, valid_links = c("logit", "log"),
  family = c("binomial", "Gamma")) {
  linktemp1 <- substitute(link1, parent.frame())
  linktemp2 <- substitute(link2, parent.frame())
  if (!is.character(linktemp1)) linktemp1 <- deparse(linktemp1)
  if (!is.character(linktemp2)) linktemp2 <- deparse(linktemp2)

  okLinks <- valid_links[1]
  if (linktemp1 %in% okLinks)
    stats1 <- stats::make.link(linktemp1)
  else if (is.character(link1))
    stats1 <- stats::make.link(link1)

  okLinks <- valid_links[2]
  if (linktemp2 %in% okLinks)
    stats2 <- stats::make.link(linktemp2)
  else if (is.character(link1))
    stats2 <- stats::make.link(link2)

  list(
    family = family,
    link = c(linktemp1, linktemp2),
    linkfun = list(stats1$linkfun, stats2$linkfun),
    linkinv = list(stats1$linkinv, stats2$linkinv),
    delta = TRUE
  )
}

#' @param link1 Link for first part of delta/hurdle model.
#' @param link2 Link for second part of delta/hurdle model.
#' @export
#' @examples
#' delta_gamma()
#' @rdname families
delta_gamma <- function(link1 = "logit", link2 = "log") {
  .make_delta(link1, link2, family = c("binomial", "Gamma"))
}

#' @export
#' @examples
#' delta_lognormal()
#' @rdname families
delta_lognormal <- function(link1 = "logit", link2 = "log") {
  .make_delta(link1, link2, family = c("binomial", "lognormal"))
}

#' @export
#' @examples
#' delta_truncated_nbinom2()
#' @rdname families
delta_truncated_nbinom2 <- function(link1 = "logit", link2 = "log") {
  .make_delta(link1, link2, family = c("binomial", "truncated_nbinom2"))
}

#' @export
#' @examples
#' delta_truncated_nbinom1()
#' @rdname families
delta_truncated_nbinom1 <- function(link1 = "logit", link2 = "log") {
  .make_delta(link1, link2, family = c("binomial", "truncated_nbinom1"))
}

#' @examples
#' delta_gamma2()
#' @rdname families
#' @details `delta_poisson_link_gamma()` is the Poisson-link (complementary
#'   log-log) delta model (Thorson 2018).
#' @references
#' Thorson, J. T. (2018). Three problems with the conventional delta-model for
#' biomass sampling data, and a computationally efficient alternative. Canadian
#' Journal of Fisheries and Aquatic Sciences, 75(9), 1369-1382.
#' \doi{10.1139/cjfas-2017-0266}
#' @export
delta_poisson_link_gamma <- function(link1 = "log", link2 = "log") {
  out <- .make_delta(link1, link2, family = c("binomial", "Gamma"),
    valid_links = c("log", "log"))
  out$type <- "poisson_link_delta"
  out
}

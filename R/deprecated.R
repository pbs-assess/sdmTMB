#' Extract MCMC samples from a model fit with \pkg{tmbstan}.
#'
#' @param object **Deprecated** See the \pkg{sdmTMBextra} package.
#'
#' @return
#' Deprecated See \pkg{sdmTMBextra} package.
#'
#' @export
extract_mcmc <- function(object = deprecated()) {
  deprecate_stop("0.2.2", "extract_mcmc()", details = xtra_msg)
}

xtra_msg <- "Load the sdmTMBextra package to use this function: https://github.com/pbs-assess/sdmTMBextra"

#' DHARMa residuals
#'
#' @param simulated_response **Deprecated** See the \pkg{sdmTMBextra} package.
#' @param object **Deprecated** See the \pkg{sdmTMBextra} package.
#' @param plot **Deprecated** See the \pkg{sdmTMBextra} package.
#' @param ... **Deprecated** See the \pkg{sdmTMBextra} package.
#'
#' @return
#' Deprecated See \pkg{sdmTMBextra} package.
#' @export

dharma_residuals <- function(simulated_response = deprecated(), object, plot = TRUE, ...) {
  deprecate_stop("0.2.2", "dharma_residuals()", details = xtra_msg)
}


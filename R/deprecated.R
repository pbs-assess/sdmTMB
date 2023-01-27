#' Extract MCMC samples from a model fit with \pkg{tmbstan}.
#'
#' Moved to the \href{https://github.com/pbs-assess/sdmTMBextra}{sdmTMBextra}
#' package
#'
#' @param object **Deprecated** See the \href{https://github.com/pbs-assess/sdmTMBextra}{sdmTMBextra} package.
#'
#' @return
#' Deprecated. See the
#' \href{https://github.com/pbs-assess/sdmTMBextra}{sdmTMBextra} package.
#'
#' @keywords internal
#' @export
#' @examples
#' \dontrun{
#' extract_mcmc()
#' }
extract_mcmc <- function(object = deprecated()) {
  deprecate_stop("0.2.2", "extract_mcmc()", details = "Load the sdmTMBextra package to use this function: https://github.com/pbs-assess/sdmTMBextra")
}

#' DHARMa residuals
#'
#' Moved to the \href{https://github.com/pbs-assess/sdmTMBextra}{sdmTMBextra}
#' package
#'
#' @param simulated_response **Deprecated** See the \pkg{sdmTMBextra} package.
#' @param object **Deprecated** See the \pkg{sdmTMBextra} package.
#' @param plot **Deprecated** See the \pkg{sdmTMBextra} package.
#' @param ... **Deprecated** See the \pkg{sdmTMBextra} package.
#'
#' @return
#' Deprecated. See the
#' \href{https://github.com/pbs-assess/sdmTMBextra}{sdmTMBextra} package.
#' @keywords internal
#' @export
#' @examples
#' \dontrun{
#' dharma_residuals()
#' }

dharma_residuals <- function(simulated_response = deprecated(), object, plot = TRUE, ...) {
  deprecate_stop("0.2.2", "dharma_residuals()", details = "Load the sdmTMBextra package to use this function: https://github.com/pbs-assess/sdmTMBextra")
}


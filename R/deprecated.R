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

#' Transform a mesh object into a mesh with correlation barriers
#'
#' @param spde_obj **Deprecated** See the \pkg{sdmTMBextra} package.
#' @param barrier_sf **Deprecated** See the \pkg{sdmTMBextra} package.
#' @param range_fraction **Deprecated** See the \pkg{sdmTMBextra} package.
#' @param proj_scaling **Deprecated** See the \pkg{sdmTMBextra} package.
#' @param plot **Deprecated** See the \pkg{sdmTMBextra} package.
#'
#' @return
#' Deprecated See \pkg{sdmTMBextra} package.
#'
#' @export

add_barrier_mesh <- function(spde_obj = deprecated(), barrier_sf,
  range_fraction = 0.2,
  proj_scaling = 1, plot = FALSE) {
  deprecate_stop("0.2.2", "add_barrier_mesh()", details = xtra_msg)
}

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


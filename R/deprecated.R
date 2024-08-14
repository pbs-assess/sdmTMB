#' Extract MCMC samples from a model fit with \pkg{tmbstan}.
#'
#' Moved to the \href{https://github.com/pbs-assess/sdmTMBextra}{sdmTMBextra}
#' package
#'
#' @param object **Deprecated** See the \href{https://github.com/pbs-assess/sdmTMBextra}{sdmTMBextra} package. Make sure to load \pkg{sdmTMBextra} *after* \pkg{sdmTMB}.
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

#' Transform a mesh object into a mesh with correlation barriers
#'
#' Moved to the \href{https://github.com/pbs-assess/sdmTMBextra}{sdmTMBextra}
#' package. Make sure to load \pkg{sdmTMBextra} *after* \pkg{sdmTMB}.
#'
#' @param spde_obj Output from [make_mesh()].
#' @param barrier_sf An sf object with polygons defining the barriers. For
#'   example, a coastline dataset for ocean data. **Note that this object must
#'   have the same projection as the data used to generate the x and y columns
#'   in `spde_obj`.**
#' @param range_fraction The fraction of the spatial range that barrier
#'   triangles have.
#' @param proj_scaling If `spde_obj` was created with scaling of the coordinates
#'   after the projection (e.g., dividing UTMs by 1000 so the spatial range is
#'   on a reasonable scale) the x and y values in `spde_obj` are multiplied by
#'   this scaling factor before applying the projection from `barrier_sf`.
#' @param plot Logical.
#' @export
#' @keywords internal
#' @return
#' Deprecated. See the
#' \href{https://github.com/pbs-assess/sdmTMBextra}{sdmTMBextra} package.
#' @examples
#' \dontrun{
#' add_barrier_mesh()
#' }
add_barrier_mesh <- function(spde_obj = deprecated(), barrier_sf = deprecated(), range_fraction = 0.2,
  proj_scaling = 1, plot = FALSE) {
  deprecate_stop("0.3.0", "add_barrier_mesh()", details = "Load the sdmTMBextra package to use this function: https://github.com/pbs-assess/sdmTMBextra \n We may move this back to the main package once the functionality joins fmesher. For now, this function requires INLA.")
}

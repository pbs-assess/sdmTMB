#' Prior distributions
#'
#' @description
#' Optional priors/penalties on model parameters.
#' This tesults in penalized likelihood within TMB
#' or can be used as priors if the model is passed
#' to \pkg{tmbstan} (see the example in [extract_mcmc()]). Note that as of 2021-07-07 Jacobian adjustments
#' **are not** made when passing to Stan for MCMC sampling.
#'
#' @description
#' `normal()` and `halfnormal()` define normal and half-normal
#' priors that, at this point, must have a location (mean)
#' parameter of 0.
#'
#' @param location Location parameter.
#' @param scale Scale parameter (SD, not variance, for normal/half-normal).
#' @param matern_range
#' @param matern_sd
#' @param prob
#' @export
#' @rdname priors
#' @examples
#' normal(0, 1)
normal <- function(location = 0, scale = 1) {
  assert_that(all(location[!is.na(location)] == 0))
  assert_that(all(scale[!is.na(scale)] > 0))
  assert_that(length(location) == length(scale))
  assert_that(sum(is.na(location)) == sum(is.na(scale)))
  c(location, scale)
}

#' @export
#' @rdname priors
#' @examples
#' halfnormal(0, 1)
halfnormal <- function(location = 0, scale = 1) {
  normal(location, scale)
}

#' @param range_gt A value one expects the spatial or spatiotemporal range is **g**reater **t**han with `1 - range_prob` probability.
#' @param sigma_lt A value one expects the spatial or spatiotemporal marginal standard deviation (`sigma_O` or `sigma_E` internally) is **l**ess **t**han with `1 - sigma_prob` probability.
#' @param range_prob Probability. See description for `range_gt`.
#' @param sigma_prob Probability. See description for `sigma_lt`.
#' @export
#' @rdname priors
#'
#' @description
#' `pc_matern()` is the Penalized Complexity prior for the Matern
#' covariance function.
#'
#' @examples
#' pc_matern(range_gt = 5, sigma_lt = 1)
#'
#' d <- subset(pcod, year > 2011)
#' pcod_spde <- make_mesh(d, c("X", "Y"), cutoff = 30)
#'
#' \dontrun{
#' # - no priors on population-level effects (`b`)
#' # - halfnormal(0, 10) prior on dispersion parameter `phi`
#' # - Matern PC priors on spatial `matern_s` and spatiotemporal
#' #   `matern_st` random field parameters
#' m <- sdmTMB(density ~ s(depth, k = 3),
#'   data = d, spde = pcod_spde, family = tweedie(),
#'   share_range = FALSE, time = "year",
#'   priors = list(
#'     phi = halfnormal(0, 10),
#'     matern_s = pc_matern(range_gt = 5, sigma_lt = 1),
#'     matern_st = pc_matern(range_gt = 5, sigma_lt = 1)
#'   )
#' )
#'
#' # - no prior on intercept
#' # - normal(0, 1) prior on depth coefficient
#' # - no prior on the dispersion parameter `phi`
#' # - Matern PC prior
#' m <- sdmTMB(density ~ depth_scaled,
#'   data = d, spde = pcod_spde, family = tweedie(),
#'   spatial_only = TRUE,
#'   priors = list(
#'     b = normal(c(NA, 0), c(NA, 1)),
#'     matern_s = pc_matern(range_gt = 5, sigma_lt = 1)
#'   )
#' )
#' }
pc_matern <- function(range_gt, sigma_lt, range_prob = 0.05, sigma_prob = 0.05) {
  c(range, range_lt_prob, sigma, sigma_gt_prob)
}

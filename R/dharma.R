#' DHARMa residuals
#'
#' Plot (and possibly return) DHARMa residuals. This is a wrapper function
#' around [DHARMa::createDHARMa()] to facilitate its use with [sdmTMB()] models.
#' **Note:** It is recommended to set `type = "mle-mvn"` in
#' [sdmTMB::simulate.sdmTMB()] for the resulting residuals to have the
#' expected distribution. This is *not* the default.
#'
#' @param simulated_response Output from [simulate.sdmTMB()]. It is recommended
#'   to set `type = "mle-mvn"` in the call to [simulate.sdmTMB()] for the
#'   residuals to have the expected distribution.
#' @param object Output from [sdmTMB()].
#' @param return_DHARMa Logical.
#' @param plot Logical.
#' @param expected_distribution Experimental: expected distribution for
#'   comparison: uniform(0, 1) or normal(0, 1). Traditional \pkg{DHARMa}
#'   residuals are uniform. If `"normal"`, a `pnorm()` transformation is applied.
#'   First, any simulated quantiles of 0 (no simulations were smaller than the
#'   observation) are set to an arbitrary value of `1/(n*10)` where `n` is the
#'   number of simulated replicated. Any simulated quantiles of 1 (no
#'   simulations were larger than the observation) are set to an arbitrary value
#'   of `1 - 1/(n*10)`. These points are shown with crosses overlaid.
#' @param ... Other arguments to pass to [DHARMa::createDHARMa()].
#'
#' @details
#' Advantages to these residuals over the ones from the [residuals.sdmTMB()]
#' method are (1) they work with delta/hurdle models for the combined
#' predictions, not the just the two parts separately, (2) they should work for
#' all families, not the just the families where we have worked out the
#' analytical quantile function, and (3) they can be used with the various
#' diagnostic tools and plots from the \pkg{DHARMa} package.
#'
#' Disadvantages are (1) they are slower to calculate since one must first
#' simulate from the model, (2) the stability of the distribution of the
#' residuals depends on having a sufficient number of simulation draws, (3)
#' uniformly distributed residuals put less emphasis on the tails visually
#' (which or may not be desired).
#'
#' Note that \pkg{DHARMa} returns residuals that are uniform(0, 1) if the data
#' are consistent with the model whereas any randomized quantile residuals from
#' [residuals.sdmTMB()] are expected to be normal(0, 1). An experimental option
#' `expected_distribution` is included to transform the distributions to
#' a normal(0, 1) expectation.
#'
#' @return
#' A data frame of observed and expected values is invisibly returned,
#' so you can set `plot = FALSE` and assign the output to an object if you wish
#' to plot the residuals yourself. See the examples.
#'
#' If `return_DHARMa = TRUE`, the object from `DHARMa::createDHARMa()`
#' is returned and any subsequent \pkg{DHARMa} functions can be applied.
#'
#' @export
#'
#' @importFrom assertthat assert_that
#' @importFrom cli cli_abort
#'
#' @seealso [simulate.sdmTMB()], [residuals.sdmTMB()]
#'
#' @examplesIf requireNamespace("DHARMa", quietly = TRUE)
#' # Try Tweedie family:
#' fit <- sdmTMB(density ~ as.factor(year) + s(depth, k = 3),
#'   data = pcod_2011, mesh = pcod_mesh_2011,
#'   family = tweedie(link = "log"), spatial = "on")
#'
#' # The `simulated_response` argument is first so the output from
#' # simulate() can be piped to dharma_residuals().
#'
#' # We will work with 100 simulations for fast examples, but you'll
#' # likely want to work with more than this (enough that the results
#' # are stable from run to run).
#'
#' # not great:
#' set.seed(123)
#' simulate(fit, nsim = 100, type = "mle-mvn") |>
#'   dharma_residuals(fit)
#'
#' \donttest{
#' # delta-lognormal looks better:
#' set.seed(123)
#' fit_dl <- update(fit, family = delta_lognormal())
#' simulate(fit_dl, nsim = 100, type = "mle-mvn") |>
#'   dharma_residuals(fit)
#'
#' # or skip the pipe:
#' set.seed(123)
#' s <- simulate(fit_dl, nsim = 100, type = "mle-mvn")
#' # and manually plot it:
#' r <- dharma_residuals(s, fit_dl, plot = FALSE)
#' head(r)
#' plot(r$expected, r$observed)
#' abline(0, 1)
#'
#' # return the DHARMa object and work with the DHARMa methods
#' ret <- simulate(fit_dl, nsim = 100, type = "mle-mvn") |>
#'   dharma_residuals(fit, return_DHARMa = TRUE)
#' plot(ret)
#'
#' # try normal(0, 1) residuals:
#' s <- simulate(fit_dl, nsim = 100, type = "mle-mvn")
#' dharma_residuals(s, fit, expected_distribution = "normal")
#' # note the points in the top right corner that had Inf quantiles
#' # because of pnorm(1)
#'
#' # work with the residuals themselves:
#' r <- dharma_residuals(s, fit, return_DHARMa = TRUE)
#' plot(fitted(fit), r$scaledResiduals)
#' }

dharma_residuals <- function(simulated_response, object, plot = TRUE,
  return_DHARMa = FALSE,
  expected_distribution = c("uniform", "normal"), ...) {
  if (!requireNamespace("DHARMa", quietly = TRUE)) {
    cli_abort("DHARMa must be installed to use this function.")
  }
  assert_that(inherits(object, "sdmTMB"))
  assert_that(is.logical(plot))
  assert_that(is.matrix(simulated_response))
  assert_that(nrow(simulated_response) == nrow(object$response))
  expected_distribution <- match.arg(expected_distribution)
  if (attr(simulated_response, "type") != "mle-mvn" && !is.null(object$tmb_random)) {
    cli_warn("It is recommended to use `simulate.sdmTMB(fit, type = 'mle-mvn')` if simulating for DHARMa residuals. See the description in ?residuals.sdmTMB under the types of residuals section.")
  }
  if (isTRUE(object$family$delta)) {
    y <- ifelse(!is.na(object$response[,2]),
      object$response[,2], object$response[,1])
  } else {
    y <- object$response[,1]
  }
  y <- as.numeric(y)
  fitted <- fitted(object)
  res <- DHARMa::createDHARMa(
    simulatedResponse = simulated_response,
    observedResponse = y,
    fittedPredictedResponse = fitted,
    ...
  )
  if (return_DHARMa) return(res)
  u <- res$scaledResiduals
  if (expected_distribution == "uniform") {
    n <- length(u)
    m <- seq_len(n) / (n + 1)
    z <- stats::qqplot(m, u, plot.it = FALSE)
    if (plot) {
      DHARMa::plotQQunif(
        res,
        testUniformity = FALSE,
        testOutliers = FALSE, testDispersion = FALSE
      )
    }
    return(invisible(data.frame(observed = z$y, expected = z$x)))
  } else { # normal
    .u <- u
    nsim <- ncol(simulated_response) * 10
    .u[.u == 1] <- 1 - (1 / nsim)
    .u[.u == 0] <- 1 / nsim
    resid_n01 <- stats::qnorm(.u)
    z <- stats::qqnorm(resid_n01, plot.it = plot);graphics::abline(0, 1)
    zo <- u %in% c(0, 1)
    if (plot) graphics::points(z$x[zo], z$y[zo], pch = 4)
    return(invisible(data.frame(observed = z$y, expected = z$x)))
  }
}


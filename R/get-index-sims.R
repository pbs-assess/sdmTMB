#' Calculate population index via simulation from joint precision matrix
#'
#' @param obj [predict.sdmTMB()] output with `sims > 0`.
#' @param level The confidence level.
#' @param return_sims Logical. Return simulation draws? The default is a
#'   quantile summary of those draws.
#' @param est_function Function to summarize the estimate (the expected value).
#'   `mean()` would be an alternative to `median()`.
#' @param agg_function Function to aggregate samples within each time slice.
#'
#' @export
#' @examples
#' pcod_spde <- make_mesh(pcod, c("X", "Y"), cutoff = 30)
#' m <- sdmTMB(density ~ 0 + as.factor(year) + depth_scaled + depth_scaled2,
#'   data = pcod, spde = pcod_spde, family = tweedie(link = "log"),
#'   time = "year"
#' )
#' p <- predict(m, newdata = qcs_grid, sims = 100)
#' x <- get_index_sims(p)
#' library(ggplot2)
#' ggplot(x, aes(year, est, ymin = lwr, ymax = upr)) +
#'   geom_line() +
#'   geom_ribbon(alpha = 0.4)
#' x_sims <- get_index_sims(p, return_sims = TRUE)
#' ggplot(x_sims, aes(as.factor(year), .value)) +
#'   geom_violin()
get_index_sims <- function(obj,
                           level = 0.95,
                           return_sims = FALSE,
                           est_function = stats::median,
                           agg_function = function(x) sum(exp(x))) {
  assert_that(is.matrix(obj), !is.null(attr(obj, "time")),
    msg = paste0("`obj` should be matrix output from `predict.sdmTMB()` ",
      "with the `sims > 0`."))
  assert_that(is.logical(return_sims))
  assert_that(is.function(est_function))
  assert_that(is.function(agg_function))
  assert_that(level > 0 && level < 1)

  .t <- as.numeric(rownames(obj))
  yrs <- sort(unique(.t))
  yr_indexes <- lapply(yrs, function(x) which(.t %in% x))
  out1 <- lapply(yr_indexes, function(x) {
    apply(obj[x, , drop = FALSE], 2L, agg_function)
  })
  if (return_sims) {
    out2 <- lapply(seq_along(out1), function(i) {
      ret <- data.frame(
        .time = yrs[i], .value = out1[[i]],
        .iteration = seq_along(out1[[i]])
      )
      stats::setNames(ret, c(attr(obj, "time"), ".value", ".iteration"))
    })
    return(do.call("rbind", out2))
  } else {
    out <- lapply(out1, function(x) {
      data.frame(
        est = est_function(x),
        lwr = stats::quantile(x, probs = (1 - level) / 2),
        upr = stats::quantile(x, probs = 1 - (1 - level) / 2),
        se = stats::sd(log(x)),
        log_est = mean(log(x))
      )
    })
    out <- do.call("rbind", out)
    out[[attr(obj, "time")]] <- yrs
    out <- out[, c(
      attr(obj, "time"), "est", "lwr", "upr", "log_est", "se"
    ), drop = FALSE]
    return(`row.names<-`(out, NULL))
  }
}

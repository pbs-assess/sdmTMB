#' Calculate population index via simulation from joint precision matrix
#'
#' @param fit_obj [sdmTMB()] output
#' @param pred_obj [predict.sdmTMB()] output with `sims > 0`
#' @param level Tail quantile
#' @param return_sims Logical. Return simulation draws (vs. quantile summary).
#' @param est_function Function to summarize expected value. `mean()` would be
#'   an alternative.
#' @param aggregate_function Function to summarize samples with each year.
#'
#' @export
#' @examples
#' pcod_spde <- make_mesh(pcod, c("X", "Y"), cutoff = 30)
#' m <- sdmTMB(density ~ 0 + as.factor(year) + depth_scaled + depth_scaled2,
#'   data = pcod, spde = pcod_spde, family = tweedie(link = "log"),
#'   time = "year")
#' p <- predict(m, newdata = qcs_grid, sims = 100)
#' x <- get_index_sims(m, p)
#' library(ggplot2)
#' ggplot(x, aes(year, est, ymin = lwr, ymax = upr)) +
#'   geom_line() + geom_ribbon(alpha = 0.4)
#' x <- get_index_sims(m, p, return_sims = TRUE)
#' ggplot(x, aes(as.factor(year), .value)) + geom_violin()

get_index_sims <- function(fit_obj, pred_obj, level = 0.95,
  return_sims = FALSE, est_function = stats::median,
  aggregate_function = function(x) sum(exp(x))) {
  assert_that(is.logical(return_sims))
  assert_that(is.function(est_function))
  assert_that(is.function(aggregate_function))
  assert_that(level > 0 && level < 1)
  assert_that(class(fit_obj) == "sdmTMB")
  assert_that(is.matrix(pred_obj))

  .t <- fit_obj$data[[fit_obj$time]]
  yrs <- sort(unique(.t))
  yr_indexes <- lapply(yrs, function(x) which(.t %in% x))
  out1 <- lapply(yr_indexes, function(x) {
    apply(pred_obj[x, , drop = FALSE], 2, aggregate_function)
  })
  if (return_sims) {
    out2 <- lapply(seq_along(out1), function(i) {
      ret <- data.frame(
        .time = yrs[i], .value = out1[[i]],
        .iteration = seq_along(out1[[i]])
      )
      stats::setNames(ret, c(fit_obj$time, ".value", ".iteration"))
    })
    return(do.call("rbind", out2))
  } else {
    out <- lapply(out1, function(x) {
      data.frame(
        est = est_function(x),
        lwr = stats::quantile(x, probs = (1 - level) / 2),
        upr = stats::quantile(x, probs = 1 - (1 - level) / 2)
      )
    })
    out <- do.call("rbind", out)
    out[[fit_obj$time]] <- yrs
    out <- out[, c(fit_obj$time, "est", "lwr", "upr"), drop = FALSE]
    return(`row.names<-`(out, NULL))
  }
}

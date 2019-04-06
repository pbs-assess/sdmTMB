#' Extract a relative biomass or abundance index
#'
#' @param obj Output from [predict.sdmTMB()].
#' @param value_name Name of value to extract.
#' @param bias_correct Should bias correction be implemented in
#'   [TMB::sdreport()]?
#' @examples
#' # Use a small number of knots for this example to make it fast:
#' pcod_spde <- make_spde(pcod$X, pcod$Y, n_knots = 50)
#' m <- sdmTMB(
#'  data = pcod,
#'  formula = density ~ 0 + as.factor(year) + depth_scaled + depth_scaled2,
#'  time = "year", spde = pcod_spde, family = tweedie(link = "log")
#' )
#' predictions <- predict(m, newdata = qcs_grid)
#' ind <- get_index(predictions, bias_correct = FALSE) # not bias correcting for speed
#'
#' library(ggplot2)
#' ggplot(ind, aes(year, est)) + geom_line() +
#'   geom_ribbon(aes(ymin = lwr, ymax = upr), alpha = 0.4)
#'
#' @export
get_index <- function(obj, value_name = "log_total", bias_correct = FALSE)  {
  sr <- TMB::sdreport(obj$obj, bias.correct = bias_correct)
  ssr <- summary(sr, "report")
  log_total <- ssr[row.names(ssr) == value_name, , drop = FALSE]
  row.names(log_total) <- NULL
  d <- as.data.frame(log_total)
  if (bias_correct)
    d <- d[,3:2,drop=FALSE]
  names(d) <- c("log_est", "se")
  d$est <- exp(d$log_est)
  d$lwr <- exp(d$log_est + stats::qnorm(0.025) * d$se)
  d$upr <- exp(d$log_est + stats::qnorm(0.975) * d$se)
  d$year <- sort(unique(obj$data$year))
  d[,c('year', 'est', 'lwr', 'upr', 'log_est', 'se'), drop = FALSE]
}

#' @rdname get_index
#' @export
get_cog <- function(obj, bias_correct = FALSE)  {
  sr <- TMB::sdreport(obj$obj, bias.correct = bias_correct)
  ssr <- summary(sr, "report")
  cog <- ssr[row.names(ssr) %in% c("cog_x", "cog_y"), , drop = FALSE]
  row.names(cog) <- NULL
  d <- as.data.frame(cog)
  if (bias_correct)
    d <- d[,3:2,drop=FALSE]
  names(d) <- c("est", "se")
  d$est <- d$est
  d$lwr <- d$est + stats::qnorm(0.025) * d$se
  d$upr <- d$est + stats::qnorm(0.975) * d$se
  d$year <- sort(unique(obj$data$year))
  d$coord <- c(rep("X", each = nrow(d)/2), rep("Y", each = nrow(d)/2))
  d[,c('coord', 'year', 'lwr', 'upr', 'est', 'se'), drop = FALSE]
}


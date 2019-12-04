#' Extract a relative biomass/abundance index or a center of gravity
#'
#' @param obj Output from [predict.sdmTMB()] with `return_tmb_object = TRUE`.
#' @param bias_correct Should bias correction be implemented [TMB::sdreport()]?
#' @param level The confidence level.
#'
#' @examples
#' # Use a small number of knots for this example to make it fast:
#' pcod_spde <- make_spde(pcod$X, pcod$Y, n_knots = 50)
#' m <- sdmTMB(
#'  data = pcod,
#'  formula = density ~ 0 + as.factor(year),
#'  time = "year", spde = pcod_spde, family = tweedie(link = "log")
#' )
#' predictions <- predict(m, newdata = qcs_grid, return_tmb_object = TRUE)
#' ind <- get_index(predictions, bias_correct = FALSE) # not bias correcting for speed
#'
#' library(ggplot2)
#' ggplot(ind, aes(year, est)) + geom_line() +
#'   geom_ribbon(aes(ymin = lwr, ymax = upr), alpha = 0.4)
#'
#' cog <- get_cog(predictions)
#' cog
#'
#' @export
get_index <- function(obj, bias_correct = FALSE, level = 0.95)  {
  d <- get_generic(obj, value_name = "log_total",
    bias_correct = bias_correct, level = level, trans = exp)
  names(d)[names(d) == "trans_est"] <- "log_est"
  d
}

#' @rdname get_index
#' @export
get_cog <- function(obj, bias_correct = FALSE, level = 0.95)  {
  d <- get_generic(obj, value_name = c("cog_x", "cog_y"),
    bias_correct = bias_correct, level = level, trans = I)
  d <- d[, names(d) != "trans_est", drop = FALSE]
  d$coord <- c(rep("X", each = nrow(d)/2), rep("Y", each = nrow(d)/2))
  d
}

get_generic <- function(obj, value_name, bias_correct = FALSE, level = 0.95,
  trans = I) {

  test <- suppressWarnings(tryCatch(obj$obj$report(), error = function(e) NA))
  if (all(is.na(test)))
    stop("It looks like the model was built with an older version of sdmTMB. ",
      "Please update the model with ",
      "`your_model <- sdmTMB:::update_model(your_model)` ",
      "first before running this function.", call. = FALSE)

  sr <- TMB::sdreport(obj$obj, bias.correct = bias_correct)
  conv <- get_convergence_diagnostics(sr)
  ssr <- summary(sr, "report")
  log_total <- ssr[row.names(ssr) %in% value_name, , drop = FALSE]
  row.names(log_total) <- NULL
  d <- as.data.frame(log_total)
  if (bias_correct)
    d <- d[,3:2,drop=FALSE]
  time_name <- obj$fit_obj$time
  names(d) <- c("trans_est", "se")
  d$est <- as.numeric(trans(d$trans_est))
  d$lwr <- as.numeric(trans(d$trans_est + stats::qnorm((1-level)/2) * d$se))
  d$upr <- as.numeric(trans(d$trans_est + stats::qnorm(1-(1-level)/2) * d$se))
  d[[time_name]] <- sort(unique(obj$data[[time_name]]))
  d$max_gradient <- max(conv$final_grads)
  d$bad_eig <- conv$bad_eig
  d[,c(time_name, 'est', 'lwr', 'upr', 'trans_est', 'se', 'max_gradient', 'bad_eig'),
    drop = FALSE]
}

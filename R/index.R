#' Extract a relative biomass/abundance index or a center of gravity
#'
#' @param obj Output from [predict.sdmTMB()] with `return_tmb_object = TRUE`.
#' @param bias_correct Should bias correction be implemented [TMB::sdreport()]?
#' @param level The confidence level.
#' @param ... Passed to [TMB::sdreport()].
#'
#' @seealso [get_index_sims()]
#'
#' @examples
#' \donttest{
#' if (inla_installed()) {
#' # Use a small number of knots for this example to make it fast:
#' pcod_spde <- make_mesh(pcod, c("X", "Y"), n_knots = 60, type = "kmeans")
#' m <- sdmTMB(
#'  data = pcod,
#'  formula = density ~ 0 + as.factor(year),
#'  time = "year", mesh = pcod_spde, family = tweedie(link = "log")
#' )
#' # Note `return_tmb_object = TRUE` and the prediction grid:
#' predictions <- predict(m, newdata = qcs_grid, return_tmb_object = TRUE)
#' ind <- get_index(predictions)
#'
#' if (require("ggplot2", quietly = TRUE)) {
#' ggplot(ind, aes(year, est)) + geom_line() +
#'   geom_ribbon(aes(ymin = lwr, ymax = upr), alpha = 0.4)
#' }
#'
#' cog <- get_cog(predictions)
#' cog
#' }
#' }
#'
#' @export
get_index <- function(obj, bias_correct = FALSE, level = 0.95, ...)  {
  d <- get_generic(obj, value_name = "link_total",
    bias_correct = bias_correct, level = level, trans = exp, ...)
  names(d)[names(d) == "trans_est"] <- "log_est"
  d
}

#' @rdname get_index
#' @param format Long or wide.
#' @export
get_cog <- function(obj, bias_correct = FALSE, level = 0.95, format = c("long", "wide"), ...)  {
  d <- get_generic(obj, value_name = c("cog_x", "cog_y"),
    bias_correct = bias_correct, level = level, trans = I, ...)
  d <- d[, names(d) != "trans_est", drop = FALSE]
  d$coord <- c(rep("X", each = nrow(d)/2), rep("Y", each = nrow(d)/2))
  format <- match.arg(format)
  if (format == "wide") {
    x <- d[d$coord == "X", c("est", "lwr", "upr", "se"),drop=FALSE]
    y <- d[d$coord == "Y", c("est", "lwr", "upr", "se"),drop=FALSE]
    names(x) <- paste0(names(x), "_", "x")
    names(y) <- paste0(names(y), "_", "y")
    d <- cbind(d[d$coord == "X", "year", drop=FALSE], cbind(x, y))
  }
  d
}

get_generic <- function(obj, value_name, bias_correct = FALSE, level = 0.95,
  trans = I, ...) {
  if (is.null(obj[["obj"]])) {
    stop("`obj` needs to be created with ",
      "`sdmTMB(..., return_tmb_object = TRUE).`", call. = FALSE)
  }
  test <- suppressWarnings(tryCatch(obj$obj$report(), error = function(e) NA))
  if (all(is.na(test)))
    stop("It looks like the model was built with an older version of sdmTMB. ",
      "Please update the model with ",
      "`your_model <- sdmTMB:::update_model(your_model)` ",
      "first before running this function.", call. = FALSE)

  tmb_data <- obj$pred_tmb_data
  if (value_name[1] == "link_total")
    tmb_data$calc_index_totals <- 1L
  if (value_name[1] == "cog_x")
    tmb_data$calc_cog <- 1L

  new_obj <- TMB::MakeADFun(
    data = tmb_data,
    parameters = get_pars(obj$fit_obj),
    map = obj$fit_obj$tmb_map,
    random = obj$fit_obj$tmb_random,
    DLL = "sdmTMB",
    silent = TRUE
  )
  # need to initialize the new TMB object once?
  # new_obj$fn(obj$fit_obj$model$par)
  if ("ADreportIndex" %in% names(new_obj$env)) {
    ind <- new_obj$env$ADreportIndex()
    to_split <- as.vector(unlist(ind[value_name]))
  } else {
    to_split <- NULL
  }

  sr <- TMB::sdreport(new_obj, bias.correct = bias_correct,
    bias.correct.control = list(sd = FALSE, split = to_split, nsplit = NULL), ...)
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
  d[[time_name]] <- sort(unique(obj$fit_obj$data[[time_name]]))
  d$max_gradient <- max(conv$final_grads)
  d$bad_eig <- conv$bad_eig
  d[,c(time_name, 'est', 'lwr', 'upr', 'trans_est', 'se', 'max_gradient', 'bad_eig'),
    drop = FALSE]
}

#' Extract a relative biomass/abundance index or a center of gravity
#'
#' @param obj Output from [predict.sdmTMB()] with `return_tmb_object = TRUE`.
#' @param bias_correct Should bias correction be implemented [TMB::sdreport()]?
#' @param level The confidence level.
#' @param area Grid cell area. A vector of length `newdata` from
#'   [predict.sdmTMB()] or a value of length 1, which will be repeated
#'   internally to match.
#' @param silent Silent?
#' @param ... Passed to [TMB::sdreport()].
#'
#' @seealso [get_index_sims()]
#' @return
#' For `get_index()`:
#' A data frame with a columns for time, estimate, lower and upper
#' confidence intervals, log estimate, and standard error of the log estimate.
#'
#' For `get_cog()`:
#' A data frame with a columns for time, estimate (center of gravity in x and y
#' coordinates), lower and upper confidence intervals, and standard error of
#' center of gravity coordinates.
#'
#' @references
#'
#' Geostatistical random-field model-based indices of abundance
#' (along with many newer papers):
#'
#' Shelton, A.O., Thorson, J.T., Ward, E.J., and Feist, B.E. 2014. Spatial
#' semiparametric models improve estimates of species abundance and
#' distribution. Canadian Journal of Fisheries and Aquatic Sciences 71(11):
#' 1655--1666. \doi{10.1139/cjfas-2013-0508}
#'
#' Thorson, J.T., Shelton, A.O., Ward, E.J., and Skaug, H.J. 2015.
#' Geostatistical delta-generalized linear mixed models improve precision for
#' estimated abundance indices for West Coast groundfishes. ICES J. Mar. Sci.
#' 72(5): 1297–1310. \doi{10.1093/icesjms/fsu243}
#'
#' Geostatistical model-based centre of gravity:
#'
#' Thorson, J.T., Pinsky, M.L., and Ward, E.J. 2016. Model-based inference for
#' estimating shifts in species distribution, area occupied and centre of
#' gravity. Methods Ecol Evol 7(8): 990–1002. \doi{10.1111/2041-210X.12567}
#'
#' Bias correction:
#'
#' Thorson, J.T., and Kristensen, K. 2016. Implementing a generic method for
#' bias correction in statistical models using random effects, with spatial and
#' population dynamics examples. Fisheries Research 175: 66–74.
#' \doi{10.1016/j.fishres.2015.11.016}
#'
#' @examplesIf inla_installed()
#' \donttest{
#' # Use a small number of knots for this example to make it fast:
#' pcod_spde <- make_mesh(pcod, c("X", "Y"), n_knots = 60, type = "kmeans")
#' m <- sdmTMB(
#'  data = pcod,
#'  formula = density ~ 0 + as.factor(year),
#'  time = "year", mesh = pcod_spde, family = tweedie(link = "log")
#' )
#'
#' # make prediction grid:
#' nd <- replicate_df(qcs_grid, "year", unique(pcod$year))
#'
#' # Note `return_tmb_object = TRUE` and the prediction grid:
#' predictions <- predict(m, newdata = nd, return_tmb_object = TRUE)
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
#'
#' @export
get_index <- function(obj, bias_correct = FALSE, level = 0.95, area = 1, silent = TRUE, ...)  {
  d <- get_generic(obj, value_name = "link_total",
    bias_correct = bias_correct, level = level, trans = exp, area = area, ...)
  names(d)[names(d) == "trans_est"] <- "log_est"
  d
}

#' @rdname get_index
#' @param format Long or wide.
#' @export
get_cog <- function(obj, bias_correct = FALSE, level = 0.95, format = c("long", "wide"), area = 1, silent = TRUE, ...)  {
  if (bias_correct) cli_abort("Bias correction with get_cog() is temporarily disabled.")
  d <- get_generic(obj, value_name = c("cog_x", "cog_y"),
    bias_correct = bias_correct, level = level, trans = I, area = area, ...)
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
  trans = I, area = 1, silent = TRUE, ...) {

  if ((!isTRUE(obj$do_index) && value_name[1] == "link_total") || value_name[1] == "cog_x") {
    if (is.null(obj[["obj"]])) {
      cli_abort(paste0("`obj` needs to be created with ",
        "`predict(..., return_tmb_object = TRUE).`"))
    }
    test <- suppressWarnings(tryCatch(obj$obj$report(obj$obj$env$last.par),
      error = function(e) NA))
    if (all(is.na(test)))
      cli_abort(c("It looks like the model was built with an older version of sdmTMB. ",
        "Please refit with the current version."))

    if (bias_correct && obj$fit_obj$control$parallel > 1) {
      cli_warn("Bias correction can be slower with multiple cores; using 1 core.")
      obj$fit_obj$control$parallel <- 1L
    }

    # FIXME parallel setup here?

    predicted_time <- sort(unique(obj$data[[obj$fit_obj$time]]))
    fitted_time <- sort(unique(obj$fit_obj$data[[obj$fit_obj$time]]))
    if (!all(fitted_time %in% predicted_time)) {
      cli_abort(paste0("Some of the fitted time elements were not predicted ",
        "on with `predict.sdmTMB()`. Please include all time elements."))
    }

    if (length(area) != nrow(obj$pred_tmb_data$proj_X_ij[[1]]) && length(area) != 1L) {
      cli_abort("`area` should be of the same length as `nrow(newdata)` or of length 1.")
    }
    if (length(area) == 1L)
      area <- rep(area, nrow(obj$pred_tmb_data$proj_X_ij[[1]]))

    tmb_data <- obj$pred_tmb_data
    tmb_data$area_i <- area
    if (value_name[1] == "link_total")
      tmb_data$calc_index_totals <- 1L
    if (value_name[1] == "cog_x")
      tmb_data$calc_cog <- 1L

    pars <- get_pars(obj$fit_obj)

    eps_name <- "eps_index" # FIXME break out into function; add for COG?
    pars[[eps_name]] <- numeric(0)

    new_obj <- TMB::MakeADFun(
      data = tmb_data,
      parameters = pars,
      map = obj$fit_obj$tmb_map,
      random = obj$fit_obj$tmb_random,
      DLL = "sdmTMB",
      silent = silent
    )

    old_par <- obj$fit_obj$model$par
    new_obj$fn(old_par) # (sometimes) need to initialize the new TMB object once!

    sr <- TMB::sdreport(new_obj, bias.correct = FALSE, ...)
  } else {
    sr <- obj$sd_report # already done in sdmTMB(do_index = TRUE)
    pars <- get_pars(obj)
    tmb_data <- obj$tmb_data
    obj <- list(fit_obj = obj) # to match regular format
    eps_name <- "eps_index"
  }
  sr_est <- as.list(sr, "Estimate", report = TRUE)

  if (bias_correct && value_name[1] == "link_total") {
    # extract and modify parameters
    pars[[eps_name]] <- rep(0, length(sr_est$total))
    new_values <- rep(0, length(sr_est$total))
    names(new_values) <- rep(eps_name, length(new_values))
    fixed <- c(obj$fit_obj$model$par, new_values)
    new_obj2 <- TMB::MakeADFun(
      data = tmb_data,
      parameters = pars,
      map = obj$fit_obj$tmb_map,
      random = obj$fit_obj$tmb_random,
      DLL = "sdmTMB",
      silent = silent
    )
    gradient <- new_obj2$gr(fixed)
    corrected_vals <- gradient[names(fixed) == eps_name]
  } else {
    if (value_name[1] == "link_total")
      cli_inform(c("Bias correction is turned off.", "
        It is recommended to turn this on for final inference."))
  }

  # # need to initialize the new TMB object once?
  # # new_obj$fn(obj$fit_obj$model$par)
  # if ("ADreportIndex" %in% names(new_obj$env)) {
  #   ind <- new_obj$env$ADreportIndex()
  #   to_split <- as.vector(unlist(ind[value_name]))
  # } else {
  #   to_split <- NULL
  # }
  #
  # sr <- TMB::sdreport(new_obj, bias.correct = bias_correct,
  #   bias.correct.control = list(sd = FALSE, split = to_split, nsplit = NULL), ...)
  conv <- get_convergence_diagnostics(sr)
  ssr <- summary(sr, "report")
  log_total <- ssr[row.names(ssr) %in% value_name, , drop = FALSE]
  row.names(log_total) <- NULL
  d <- as.data.frame(log_total)
  # if (bias_correct)
  #   d <- d[,3:2,drop=FALSE]
  time_name <- obj$fit_obj$time
  names(d) <- c("trans_est", "se")
  if (bias_correct) {
    d$trans_est <- log(corrected_vals)
    d$est <- corrected_vals
  } else {
    d$est <- as.numeric(trans(d$trans_est))
  }
  d$lwr <- as.numeric(trans(d$trans_est + stats::qnorm((1-level)/2) * d$se))
  d$upr <- as.numeric(trans(d$trans_est + stats::qnorm(1-(1-level)/2) * d$se))
  d[[time_name]] <- sort(unique(obj$fit_obj$data[[time_name]]))
  # d$max_gradient <- max(conv$final_grads)
  # d$bad_eig <- conv$bad_eig
  d[,c(time_name, 'est', 'lwr', 'upr', 'trans_est', 'se'),
    drop = FALSE]
}

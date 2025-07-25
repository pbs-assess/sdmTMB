#' Extract a relative biomass/abundance index, center of gravity, or effective
#' area occupied
#'
#' @param obj Output from [predict.sdmTMB()] with `return_tmb_object = TRUE`.
#'   Alternatively, if [sdmTMB()] was called with `do_index = TRUE` or if using
#'   the helper function [get_index_split()], an object from [sdmTMB()].
#' @param bias_correct Should bias correction be implemented [TMB::sdreport()]?
#'   This is recommended to be `TRUE` for any final analyses, but one may wish
#'   to set this to `FALSE` for slightly faster calculations while experimenting
#'   with models.
#' @param level The confidence level.
#' @param area Grid cell area. A vector of length `newdata` from
#'   [predict.sdmTMB()] *or* a value of length 1 which will be repeated
#'   internally to match *or* a character value representing the column
#'   used for area weighting.
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
#' For `get_eao()`:
#' A data frame with a columns for time, estimate (effective area occupied; EAO),
#' lower and upper confidence intervals,
#' log EAO, and standard error of the log EAO estimates.
#'
#' @references
#'
#' Geostatistical model-based indices of abundance
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
#' Geostatistical model-based effective area occupied:
#'
#' Thorson, J.T., Rindorf, A., Gao, J., Hanselman, D.H., and Winker, H. 2016.
#' Density-dependent changes in effective area occupied for
#' sea-bottom-associated marine fishes. Proceedings of the Royal Society B:
#' Biological Sciences 283(1840): 20161853.
#' \doi{10.1098/rspb.2016.1853}
#'
#' Bias correction:
#'
#' Thorson, J.T., and Kristensen, K. 2016. Implementing a generic method for
#' bias correction in statistical models using random effects, with spatial and
#' population dynamics examples. Fisheries Research 175: 66–74.
#' \doi{10.1016/j.fishres.2015.11.016}
#'
#' @examplesIf ggplot2_installed()
#' \donttest{
#' library(ggplot2)
#'
#' # use a small number of knots for this example to make it fast:
#' mesh <- make_mesh(pcod, c("X", "Y"), n_knots = 60)
#'
#' # fit a spatiotemporal model:
#' m <- sdmTMB(
#'  data = pcod,
#'  formula = density ~ 0 + as.factor(year),
#'  time = "year", mesh = mesh, family = tweedie(link = "log")
#' )
#'
#' # prepare a prediction grid:
#' nd <- replicate_df(qcs_grid, "year", unique(pcod$year))
#'
#' # Note `return_tmb_object = TRUE` and the prediction grid:
#' predictions <- predict(m, newdata = nd, return_tmb_object = TRUE)
#'
#' # biomass index:
#' ind <- get_index(predictions, bias_correct = TRUE)
#' ind
#' ggplot(ind, aes(year, est)) + geom_line() +
#'   geom_ribbon(aes(ymin = lwr, ymax = upr), alpha = 0.4) +
#'   ylim(0, NA)
#'
#' # do that in 2 chunks
#' # only necessary for very large grids to save memory
#' # will be slower but save memory
#' # note the first argument is the model fit object:
#' ind <- get_index_split(m, newdata = nd, nsplit = 2, bias_correct = TRUE)
#'
#' # center of gravity:
#' cog <- get_cog(predictions, format = "wide")
#' cog
#' ggplot(cog, aes(est_x, est_y, colour = year)) +
#'   geom_point() +
#'   geom_linerange(aes(xmin = lwr_x, xmax = upr_x)) +
#'   geom_linerange(aes(ymin = lwr_y, ymax = upr_y)) +
#'   scale_colour_viridis_c()
#'
#' # effective area occupied:
#' eao <- get_eao(predictions)
#' eao
#' ggplot(eao, aes(year, est)) + geom_line() +
#'   geom_ribbon(aes(ymin = lwr, ymax = upr), alpha = 0.4) +
#'   ylim(0, NA)
#' }
#' @export
get_index <- function(obj, bias_correct = TRUE, level = 0.95, area = 1, silent = TRUE, ...)  {
  # if offset is a character vector, use the value in the dataframe
  if (is.character(area)) {
    area <- obj$data[[area]]
  }

  d <- get_generic(obj, value_name = "link_total",
    bias_correct = bias_correct, level = level, trans = exp, area = area, ...)
  names(d)[names(d) == "trans_est"] <- "log_est"
  d$type <- "index"
  d
}

chunk_time <- function(x, chunks) {
  assert_that(is.numeric(chunks))
  assert_that(chunks > 0)
  chunks <- as.integer(chunks)
  if (chunks == 1L) {
    list(x)
  } else {
    return(split(x, cut(seq_along(x), chunks, labels = FALSE)))
  }
}

#' @rdname get_index
#' @param newdata New data (e.g., a prediction grid by year) to pass to
#'   [predict.sdmTMB()] in the case of `get_index_split()`.
#' @param nsplit The number of splits to do the calculation in. For memory
#'   intensive operations (large grids and/or models), it can be helpful to
#'   do the prediction, area integration, and bias correction on subsets of
#'   time slices (e.g., years) instead of all at once. If `nsplit > 1`, this
#'   will usually be slower but with reduced memory use.
#' @param predict_args A list of arguments to pass to [predict.sdmTMB()] in the
#'   case of `get_index_split()`.
#' @export
get_index_split <- function(
    obj, newdata, bias_correct = FALSE, nsplit = 1,
    level = 0.95, area = 1, silent = FALSE, predict_args = list(), ...) {
  if (!inherits(obj, "sdmTMB")) {
    cli_abort("get_index_split() is meant to be run on a fitted object from sdmTMB() and not a prediction object as in get_index().")
  }
  assert_that(is.list(predict_args))

  predict_args[["return_tmb_object"]] <- TRUE
  predict_args[["object"]] <- obj

  times <- sort(obj$time_lu$time_from_data)
  time_chunks <- chunk_time(times, nsplit)

  if ("offset" %in% names(predict_args)) {
    if (!is.numeric(predict_args$offset))
      cli_abort("`offset` should be a numeric vector for use with `get_index_split()`")
    offset <- predict_args$offset
    predict_args$offset <- NULL
  } else {
    offset <- rep(0, nrow(newdata))
  }
  if (length(area) == 1L) area <- rep(area, nrow(newdata))

  msg <- paste0("Calculating index in ", nsplit, " chunks")
  if (!silent) cli::cli_progress_bar(msg, total = length(time_chunks))
  index_list <- list()
  for (i in seq_along(time_chunks)) {
    if (!silent) {
      cli::cli_progress_update(set = i, total = length(time_chunks), force = TRUE)
    }
    this_chunk_i <- newdata[[obj$time]] %in% time_chunks[[i]]
    nd <- newdata[this_chunk_i, , drop = FALSE]

    predict_args[["newdata"]] <- nd
    predict_args[["offset"]] <- offset[this_chunk_i]
    pred <- do.call(predict, predict_args)
    index_list[[i]] <-
      get_index(
        pred,
        bias_correct = bias_correct,
        level = level,
        area = area[this_chunk_i],
        silent = TRUE,
        ...
      )
  }
  if (!silent) cli::cli_progress_done()
  do.call(rbind, index_list)
}

#' @rdname get_index
#' @param format Long or wide.
#' @export
get_cog <- function(obj, bias_correct = FALSE, level = 0.95, format = c("long", "wide"), area = 1, silent = TRUE, ...)  {
  # if offset is a character vector, use the value in the dataframe
  if (is.character(area)) {
    area <- obj$data[[area]]
  }
  pred_time <- sort(unique(obj$data[[obj$fit_obj$time]]))
  fitted_time <- obj$fit_obj$fitted_time
  if (bias_correct && sum(!fitted_time %in% pred_time) > 0L)
    cli_abort("Please include all time elements in the prediction data frame if using bias_correct = TRUE with get_cog().")

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
    d <- cbind(d[d$coord == "X", obj$fit_obj$time, drop=FALSE], cbind(x, y))
  }
  d$type <- "cog"
  d
}

#' @rdname get_index
#' @export
get_eao <- function(obj,
  bias_correct = FALSE,
  level = 0.95,
  area = 1,
  silent = TRUE,
  ...
)  {
  # if offset is a character vector, use the value in the dataframe
  if (is.character(area)) {
    area <- obj$data[[area]]
  }

  d <- get_generic(obj, value_name = c("log_eao"),
    bias_correct = bias_correct, level = level, trans = exp, area = area, ...)
  names(d)[names(d) == "trans_est"] <- "log_est"
  d$type <- "eoa"
  d
}

get_generic <- function(obj, value_name, bias_correct = FALSE, level = 0.95,
  trans = I, area = 1, silent = TRUE, ...) {

  # if offset is a character vector, use the value in the dataframe
  if (is.character(area)) {
    area <- obj$data[[area]]
  }

  reinitialize(obj$fit_obj)

  if ((!isTRUE(obj$do_index) && value_name[1] == "link_total") || value_name[1] == "cog_x" || value_name[[1]] == "log_eao") {
    if (is.null(obj[["obj"]])) {
      cli_abort(paste0("`obj` needs to be created with ",
        "`predict(..., return_tmb_object = TRUE).`"))
    }

    if (!"report" %in% names(obj$obj)) {
      cli_abort(c("It looks like the predict function was run without `newdata` specified.",
        "Re-run the predict function with `newdata` specified.")) #276
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

    assert_that(!is.null(area))
    if (length(area) > 1L) {
      n_fakend <- if (!is.null(obj$fake_nd)) nrow(obj$fake_nd) else 0L
      area <- c(area, rep(1, n_fakend)) # pad area with any extra time
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
    if (value_name[1] == "log_eao")
      tmb_data$calc_eao <- 1L

    pars <- get_pars(obj$fit_obj)

    eps_name <- "eps_index" # FIXME break out into function; add for COG?
    pars[[eps_name]] <- numeric(0)

    new_obj <- TMB::MakeADFun(
      data = tmb_data,
      parameters = pars,
      profile = obj$fit_obj$control$profile,
      map = obj$fit_obj$tmb_map,
      random = obj$fit_obj$tmb_random,
      DLL = "sdmTMB",
      silent = silent
    )

    old_par <- obj$fit_obj$model$par
    new_obj$fn(old_par) # (sometimes) need to initialize the new TMB object once!

    bc <- FALSE ## done below
    sr <- TMB::sdreport(new_obj, bias.correct = bc, ...)
  } else {
    sr <- obj$sd_report # already done in sdmTMB(do_index = TRUE)
    pars <- get_pars(obj)
    tmb_data <- obj$tmb_data
    obj <- list(fit_obj = obj) # to match regular format
    eps_name <- "eps_index"
  }
  sr_est <- as.list(sr, "Estimate", report = TRUE)

  if (bias_correct && value_name[[1]] %in% c("link_total", "cog_x")) {
    # extract and modify parameters
    if (value_name[[1]] == "link_total") .n <- length(sr_est$total)
    if (value_name[[1]] == "cog_x") .n <- length(sr_est$cog_x) * 2 # 2 b/c x and y
    pars[[eps_name]] <- rep(0, .n)
    new_values <- rep(0, .n)
    names(new_values) <- rep(eps_name, length(new_values))
    fixed <- c(obj$fit_obj$model$par, new_values)
    new_obj2 <- TMB::MakeADFun(
      data = tmb_data,
      parameters = pars,
      map = obj$fit_obj$tmb_map,
      profile = obj$fit_obj$control$profile,
      random = obj$fit_obj$tmb_random,
      DLL = "sdmTMB",
      silent = silent,
      intern = FALSE, # tested as faster for most models
      inner.control = list(sparse = TRUE, lowrank = TRUE, trace = TRUE)
    )
    gradient <- new_obj2$gr(fixed)
    corrected_vals <- gradient[names(fixed) == eps_name]
  } else {
    if (value_name[[1]] == "link_total" || value_name[[1]] == "cog_x")
      cli_inform(c("Bias correction is turned off.", "
        It is recommended to turn this on for final inference."))
  }
  ssr <- summary(sr, "report")
  log_total <- ssr[row.names(ssr) %in% value_name, , drop = FALSE]
  row.names(log_total) <- NULL
  d <- as.data.frame(log_total)
  time_name <- obj$fit_obj$time
  names(d) <- c("trans_est", "se")
  if (bias_correct) {
    if (value_name[[1]] == "cog_x") {
      d$trans_est <- corrected_vals
    } else {
      d$trans_est <- log(corrected_vals)
    }
    d$est <- corrected_vals
  } else {
    d$est <- as.numeric(trans(d$trans_est))
  }
  d$lwr <- as.numeric(trans(d$trans_est + stats::qnorm((1-level)/2) * d$se))
  d$upr <- as.numeric(trans(d$trans_est + stats::qnorm(1-(1-level)/2) * d$se))

  if ("pred_tmb_data" %in% names(obj)) { # standard case
    ii <- sort(unique(obj$pred_tmb_data$proj_year))
  } else { # fit with do_index = TRUE
    ii <- sort(unique(obj$fit_obj$tmb_data$proj_year))
  }
  d <- d[d$est != 0, ,drop=FALSE] # these were not predicted on
  d <- d[!is.na(d$est), ,drop=FALSE] # these were not predicted on
  lu <- obj$fit_obj$time_lu
  tt <- lu$time_from_data[match(ii, lu$year_i)]
  if (nrow(d) == 0L) {
    msg <- c(
      "There were no results returned by TMB.",
      "It's possible TMB ran out of memory.",
      "You could try a computer with more RAM or see the function `get_index_split()` with `nsplit > 1`",
      "which lets you split the TMB sdreport and bias correction into chunks."
    )
    cli_abort(msg)
  }
  d[[time_name]] <- tt

  # remove padded extra time fake data:
  if (!is.null(obj$fake_nd)) {
    d <- d[!d[[obj$fit_obj$time]] %in% obj$fake_nd[[obj$fit_obj$time]], ,drop = FALSE]
  }
  if ("do_index_time_missing_from_nd" %in% names(obj$fit_obj)) {
    d <- d[!d[[obj$fit_obj$time]] %in% obj$fit_obj$do_index_time_missing_from_nd, ,drop = FALSE]
  }
  d[,c(time_name, 'est', 'lwr', 'upr', 'trans_est', 'se'), drop = FALSE]
}

set_par_value <- function(opt, par) {
  as.numeric(opt$par[par == names(opt$par)])
}

#' Optimization control options
#'
#' Any arguments to pass to [stats::nlminb()].
#'
#' @param eval.max Maximum number of evaluations of the objective function
#'   allowed.
#' @param iter.max Maximum number of iterations allowed.
#' @param ... Anything else. See the 'Control parameters' section of
#'   [stats::nlminb()].
#'
#' @export
sdmTMBcontrol <- function(eval.max = 1e4, iter.max = 1e4, ...) {
  list(eval.max = eval.max, iter.max = iter.max, ...)
}

get_convergence_diagnostics <- function(sd_report) {
  final_grads <- sd_report$gradient.fixed
  bad_eig <- FALSE
  if (!is.null(sd_report$pdHess)) {
    if (!sd_report$pdHess) {
      warning("The model may not have converged: ",
        "non-positive-definite Hessian matrix.", call. = FALSE)
    } else {
      eigval <- try(1 / eigen(sd_report$cov.fixed)$values, silent = TRUE)
      if (is(eigval, "try-error") || (min(eigval) < .Machine$double.eps * 10)) {
        warning("The model may not have converged: ",
          "extreme or very small eigen values detected.", call. = FALSE)
        bad_eig <- TRUE
      }
      if (any(final_grads > 0.01))
        warning("The model may not have converged. ",
          "Maximum final gradient: ", max(final_grads), ".", call. = FALSE)
    }
  }
  invisible(list(final_grads = final_grads, bad_eig = bad_eig))
}

make_year_i <- function(x) {
  x <- as.integer(as.factor(x))
  x - min(x)
}

check_offset <- function(formula) {
  any(grepl("^offset$",
    gsub(" ", "", unlist(strsplit(as.character(formula), "\\+")))))
}

update_model <- function(object, silent = FALSE) {
  if (!"nobs_RE" %in% names(object$tmb_data)) {
    object$tmb_data$nobs_RE <- 0L
    object$tmb_data$ln_tau_G_index <- 0L
    object$tmb_data$RE_indexes <- matrix(ncol = 0L, nrow = nrow(object$tmb_data$X_ij))
    object$tmb_data$proj_RE_indexes <- matrix(ncol = 0L, nrow = 1L)
    object$tmb_params$ln_tau_G <- 0
    object$tmb_params$RE <- rep(0, 1L)
    object$tmb_map$ln_tau_G <- factor(NA)
    object$tmb_map$RE <- factor(NA)
    object$split_formula <- list()
    object$split_formula$fixedFormula <- object$formula
  }
  if (!"barrier" %in% names(object$tmb_data)) {
    object$tmb_data$barrier_scaling <- c(1, 1)
    object$tmb_data$barrier <- 0L
    C0 <- rep(1, 2)
    C1 <- rep(1, 2)
    D0 <- Matrix::Matrix(0, 1, 1, doDiag = FALSE)
    D1 <- Matrix::Matrix(0, 1, 1, doDiag = FALSE)
    .I <- Matrix::Matrix(0, 1, 1, doDiag = FALSE)
    object$tmb_data$spde_barrier <- make_barrier_spde(object$spde)
  }
  if (!"pop_pred" %in% names(object$tmb_data)) object$tmb_data$pop_pred <- 0L
  if (!"penalties" %in% names(object$tmb_data)) object$tmb_data$penalties <- rep(1, ncol(object$tmb_data$X_ij))
  if (!"mgcv" %in% names(object)) object$mgcv <- FALSE
  object$tmb_data$weights_i <- rep(1, length(object$tmb_data$y_i))
  object$tmb_data$calc_quadratic_range <- 0L
  object$tmb_data$area_i <- rep(1, length(object$tmb_data$y_i))
  if (!"X_threshold" %in% names(object$tmb_data)) {
    object$tmb_data$X_threshold <- rep(0, nrow(object$data)) # just placeholder
    object$tmb_data$threshold_func <- 0L
    object$tmb_data$proj_X_threshold <- 0 # dummy
    object$tmb_params$b_threshold <- rep(0, 2)
    object$tmb_map$b_threshold <- factor(c(NA, NA))
  }
  if (!"df" %in% names(object$tmb_data)) object$tmb_data$df <- 3
  object$tmb_obj <- TMB::MakeADFun(
    data = object$tmb_data, parameters = object$tmb_params,
    map = object$tmb_map, random = object$tmb_random, DLL = "sdmTMB", silent = silent,
    checkParameterOrder = FALSE)
  object
}

check_and_parse_thresh_params <- function(formula, data) {
  terms <- stats::terms(formula)
  terms_labels <- attr(terms, "term.labels")
  if (any(grepl("linear_thresh", terms_labels)) && any(grepl("logistic_thresh", terms_labels))) {
    stop("Please include only a linear (`breakpt`) *or* a logistic threshold.", call. = FALSE)
  }
  if (sum(grepl("linear_thresh", terms_labels)) > 1 || sum(grepl("logistic_thresh", terms_labels)) > 1) {
    stop("Please include only a *single* threshold variable.", call. = FALSE)
  }
  threshold_parameter <- NULL
  if (any(grepl("breakpt", terms_labels))) {
    out <- parse_threshold_formula(formula, "breakpt", terms_labels)
    threshold_parameter <- out$threshold_parameter
    formula <- out$formula
    threshold_function <- "linear"
  }
  if (any(grepl("logistic", terms_labels))) {
    out <- parse_threshold_formula(formula, "logistic", terms_labels)
    threshold_parameter <- out$threshold_parameter
    formula <- out$formula
    threshold_function <- "logistic"
  }
  if(!is.null(threshold_parameter)) {
    if (length(threshold_parameter) > 1) {
      stop("`threshold_parameter` must be a single variable name.", call. = FALSE)
    }
    if (!threshold_parameter %in% names(data)) {
      stop("`threshold_parameter` is not a column in the `data` data frame.", call. = FALSE)
    }
  }

  if (is.null(threshold_parameter)) {
    X_threshold <- rep(0, nrow(data)) # just placeholder
    threshold_func <- 0L
  } else {
    X_threshold <- data[, names(data) == threshold_parameter, drop = TRUE]
    # indexed 1, 2 because 0 will tell TMB not to estimate this:
    threshold_func <- match(threshold_function, c("linear", "logistic"))
  }
  X_threshold <- as.numeric(unlist(X_threshold))
  list(formula = formula, threshold_parameter = threshold_parameter,
    threshold_func = threshold_func, X_threshold = X_threshold
  )
}

parse_threshold_formula <- function(formula, thresh_type_short = "lin_thresh",
  terms_labels) {
  which_thresh <- grep(thresh_type_short, terms_labels)
  temp <- gsub(paste0("^", thresh_type_short, "\\("), "", terms_labels[which_thresh])
  threshold_parameter <- gsub("\\)$", "", temp)
  formula <- stats::update(formula, paste("~ . -", terms_labels[which_thresh]))
  list(formula = formula, threshold_parameter = threshold_parameter)
}

expand_time <- function(df, time_slices, time_column) {
  df[["weight_sdmTMB"]] <- 1
  fake_df <- df[1L, , drop = FALSE]
  fake_df[["weight_sdmTMB"]] <- 0
  missing_years <- time_slices[!time_slices %in% df[[time_column]]]
  fake_df <- do.call("rbind", replicate(length(missing_years), fake_df, simplify = FALSE))
  fake_df[[time_column]] <- missing_years
  rbind(df, fake_df)
}

safe_deparse <- function (x, collapse = " ") {
  paste(deparse(x, 500L), collapse = collapse)
}

barnames <- function (bars) {
  vapply(bars, function(x) safe_deparse(x[[3]]), "")
}

check_valid_factor_levels <- function(x, .name = "") {
  assert_that(is.factor(x),
    msg = sprintf("Random effect group column `%s` is not a factor.", .name))
  lev <- sort(levels(x))
  uni <- sort(unique(as.character(x)))
  assert_that(identical(lev, uni),
    msg = sprintf(
      "Random effect group column `%s` has extra factor levels. Please remove them.", .name))
}

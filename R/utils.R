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
  if (!"barrier" %in% names(object$tmb_data)) {
    object$tmb_data$barrier_scaling <- c(1, 1)
    object$tmb_data$barrier <- 0L
    C0 <- rep(1, 2)
    C1 <- rep(1, 2)
    D0 <- Matrix::Matrix(0, 1, 1)
    D1 <- Matrix::Matrix(0, 1, 1)
    .I <- Matrix::Matrix(0, 1, 1)
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
  }
  if (!"df" %in% names(object$tmb_data)) object$tmb_data$df <- 3
  object$tmb_obj <- TMB::MakeADFun(
    data = object$tmb_data, parameters = object$tmb_params,
    map = object$tmb_map, random = object$tmb_random, DLL = "sdmTMB", silent = silent)
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

# f <- y ~ x + (1|g)
# data <- data.frame(x = runif(100), y = rnorm(100), g = gl(10, 10))
#
safe_deparse <- function (x, collapse = " ") {
  paste(deparse(x, 500L), collapse = collapse)
}
barnames <- function (bars) {
  vapply(bars, function(x) safe_deparse(x[[3]]), "")
}
# f <- y ~ x + (1|g) + (1|h)
# ff <- glmmTMB::splitForm(f)
# ff$fixedFormula
# ff$reTrmFormulas
# barnames(ff$reTrmFormulas)

# f <- y ~ x + (1|g) + (1|h)
# data <- data.frame(x = runif(100), y = rnorm(100),
#   g = gl(10, 10, labels = sample(LETTERS[1:10], size = 10)),
#   h = gl(5, 20, labels = sample(LETTERS[1:5], size = 5))
#   )
# ff <- glmmTMB::splitForm(f)
# ff$fixedFormula
# ff$reTrmFormulas
#
# # X_ij <- model.matrix(~ g, data)
# # mf <- model.frame(formula, data)
#
# data <- data[sample(1:nrow(data), nrow(data)),]
#
# RE_names <- barnames(ff$reTrmFormulas)
# RE_indexes <- vapply(RE_names, function(x) as.integer(data[[x]]), rep(1L, nrow(data)))
#
# g_index <- RE_indexes[,1]
# h_index <- RE_indexes[,2]
#
# nobs_RE <- apply(RE_indexes, 2, max)
#
# RE_g <- rnorm(nobs_RE[1], mean = -10)
# RE_h <- rnorm(nobs_RE[2], mean = 30)
# RE <- c(RE_g, RE_h)
#
# # ObsinRE <- matrix(g_index, nrow = nrow(data), ncol = 1L)
# # ObsinRE <- cbind(ObsinRE, h_index)
# n_RE <- ncol(ObsinRE)
# mu <- rep(0, nrow(data))
# temp <- 0
# N <- nrow(data)
# sigmas <- c(0.8, 1.2)
#
# # for (int i=0; i<Nfish_LA; i++) {
# for (i in 1:N) {
#   # for (int k=0; k<n_RE_LA; k++){
#   for (k in 1:n_RE) {
#     if (k == 1) mu[i] = mu[i] + RE[RE_indexes[i, k]]
#     if (k > 1) {
#       temp = temp + nobs_RE[k - 1]
#       mu[i] = mu[i] + RE[RE_indexes[i, k] + temp]
#     }
#     # nll <- dnorm(RE[RE_indexes[i, k], 0, sigmas[k])
#   }
#   temp = 0
# }
# mu
#
#
# identical(mu, RE_g[g_index] + RE_h[h_index])

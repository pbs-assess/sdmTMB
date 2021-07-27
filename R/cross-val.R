ll_gaussian <- function(object, withheld_y, withheld_mu) {
  .sd <- exp(object$model$par[["ln_phi"]])
  stats::dnorm(x = withheld_y, mean = withheld_mu, sd = .sd, log = TRUE)
}

ll_tweedie <- function(object, withheld_y, withheld_mu) {
  p <- stats::plogis(object$model$par[["thetaf"]]) + 1
  phi <- exp(object$model$par[["ln_phi"]])
  fishMod::dTweedie(y = withheld_y, mu = withheld_mu, p = p, phi = phi, LOG = TRUE)
}

ll_binomial <- function(object, withheld_y, withheld_mu) {
  stats::dbinom(x = withheld_y, size = 1, prob = withheld_mu, log = TRUE)
}

ll_gamma <- function(object, withheld_y, withheld_mu) {
  .shape <- exp(object$model$par[["ln_phi"]])
  stats::dgamma(x = withheld_y, shape = .shape, scale = withheld_mu / .shape, log = TRUE)
}

ll_lognormal <- function(object, withheld_y, withheld_mu) {
  .sd <- exp(object$model$par[["ln_phi"]])
  stats::dlnorm(x = withheld_y, meanlog = withheld_mu - 0.5 * (.sd)^2, sdlog = .sd, log = TRUE)
}

ll_sdmTMB <- function(object, withheld_y, withheld_mu) {

  family_func <- switch(object$family$family,
    gaussian = ll_gaussian,
    tweedie = ll_tweedie,
    binomial = ll_binomial,
    lognormal = ll_lognormal,
    Gamma = ll_gamma,
    stop(object$family$family, " not yet implemented. ",
      "Please file an issue on GitHub.",
      call. = FALSE
    )
  )
  family_func(object, withheld_y, withheld_mu)
}

#' Save log likelihoods of k-fold cross-validation for sdmTMB models
#'
#' @param formula Model formula.
#' @param data A data frame.
#' @param spde Output from [make_mesh()].
#' @param time The name of the time column. Leave as `NULL` if this is only spatial data.
#' @param k_folds Number of folds.
#' @param fold_ids Optional vector containing user fold IDs. Can also be a
#'   single string, e.g. `"fold_id"` representing the name of the variable in
#'   `data`.
#' @param use_initial_fit Fit the first fold and use those parameter values
#'   as starting values for subsequent folds? Can be faster with many folds.
#' @param ... All other arguments required to run [sdmTMB()] model with the
#'   exception of `weights`, which are used to define the folds.
#'
#' @export
#'
#' @examples
#' if (inla_installed()) {
#' spde <- make_mesh(pcod, c("X", "Y"), cutoff = 25)
#'
#' m_cv <- sdmTMB_cv(
#'   density ~ 0 + depth_scaled + depth_scaled2,
#'   data = pcod, spde = spde,
#'   family = tweedie(link = "log"), k_folds = 2
#' )
#' m_cv$fold_loglik
#' m_cv$sum_loglik
#' head(m_cv$data)
#' m_cv$models[[1]]
#' m_cv$max_gradients
#' }
sdmTMB_cv <- function(formula, data, spde, time = NULL,
                      k_folds = 10, fold_ids = NULL,
                      use_initial_fit = FALSE,
                      ...) {
  if (k_folds < 1) stop("`k_folds` must be >= 1.", call. = FALSE)

  data[["_sdm_order_"]] <- seq_len(nrow(data))

  # add column of fold_ids stratified across time steps
  if (is.null(time)) {
    time <- "_sdmTMB_time"
    data[[time]] <- 0L
  }
  if (is.null(fold_ids)) {
    dd <- lapply(split(data, data[[time]]), function(x) {
      x$cv_fold <- sample(rep(seq(1L, k_folds), nrow(x)), size = nrow(x))
      x
    })
    data <- do.call(rbind, dd)
    fold_ids <- "cv_fold"
  } else {
    # fold_ids passed in; can be numeric, or a named column in `data`
    data$cv_fold <- NA
    if (length(fold_ids) == nrow(data)) {
      data$cv_fold <- fold_ids
    }
    if (length(fold_ids) == 1L && is.character(fold_ids)) {
      if (!fold_ids %in% names(data)) {
        stop("Name of fold identifier not found in data.", call. = FALSE)
      }
      data$cv_fold <- data[[fold_ids]]
    }
    if (length(fold_ids) > 1 && length(fold_ids) < nrow(data)) {
      stop("Dimension of `fold_ids` doesn't match data and is not a named variable.")
    }
    if (length(which(is.na(data$cv_fold))) > 0) {
      stop("NAs found in `fold_ids`; please check `fold_ids` are specified correctly.")
    }
  }

  dot_args <- as.list(substitute(list(...)))[-1L]
  if ("weights" %in% names(dot_args)) {
    stop("`weights` cannot be specified within sdmTMB_cv().", call. = FALSE)
  }

  if (k_folds > 1) {
    # data in kth fold get weight of 0:
    weights <- ifelse(data$cv_fold == 1L, 1, 0)
  } else {
    weights <- rep(1, nrow(data))
  }

  if (use_initial_fit) {
    # run model on first fold to get starting values:
    fit1 <- sdmTMB(
      data = data, formula = formula, time = time, spde = spde,
      weights = weights, ...
    )
  }

  fit_func <- function(k) {
    # data in kth fold get weight of 0:
    weights <- ifelse(data$cv_fold == k, 1, 0)
    args <- c(list(
      data = data, formula = formula, time = time,
      spde = spde, weights = weights, previous_fit = if (use_initial_fit) fit1 else NULL
    ), dot_args)
    if (k == 1L && use_initial_fit) {
      object <- fit1
    } else {
      object <- do.call(sdmTMB, args)
      if (max(object$gradients) > 0.01) {
        object <- run_extra_optimization(object, nlminb_loops = 1L, newton_loops = 0L)
      }
    }

    # predict for withheld data:
    predicted <- predict(object)[weights == 0, , drop = FALSE]
    cv_data <- data[weights == 0, , drop = FALSE]
    cv_data$cv_predicted <- object$family$linkinv(predicted$est)
    response <- get_response(object$formula)
    withheld_y <- predicted[[response]]
    withheld_mu <- cv_data$cv_predicted

    # calculate log likelihood for each withheld observation:
    cv_data$cv_loglik <- ll_sdmTMB(object, withheld_y, withheld_mu)

    list(
      data = cv_data,
      model = object,
      pdHess = object$sd_report$pdHess,
      max_gradient = max(abs(object$gradients)),
      bad_eig = object$bad_eig
    )
  }

  if (requireNamespace("future.apply", quietly = TRUE)) {
    message("Running fits with `future.apply()`.\n",
      "Set a parallel `future::plan()` to use parallel processing.")
    out <- future.apply::future_lapply(seq_len(k_folds), fit_func, future.seed = TRUE)
  } else {
    message("Running fits sequentially.\n",
      "Install the future and future.apply packages and\n",
      "set a parallel `future::plan()` to use parallel processing.")
    out <- lapply(seq_len(k_folds), fit_func)
  }

  models <- lapply(out, `[[`, "model")
  data <- lapply(out, `[[`, "data")
  fold_cv_ll <- vapply(data, function(.x) sum(.x$cv_loglik), FUN.VALUE = numeric(1L))
  data <- do.call(rbind, data)
  data <- data[order(data[["_sdm_order_"]]), , drop = FALSE]
  data[["_sdm_order_"]] <- NULL
  data[["_sdmTMB_time"]] <- NULL
  row.names(data) <- NULL
  bad_eig <- vapply(out, `[[`, "bad_eig", FUN.VALUE = logical(1L))
  pdHess <- vapply(out, `[[`, "pdHess", FUN.VALUE = logical(1L))
  max_grad <- vapply(out, `[[`, "max_gradient", FUN.VALUE = numeric(1L))
  converged <- all(!bad_eig) && all(pdHess)
  list(
    data = data,
    models = models,
    fold_loglik = fold_cv_ll,
    sum_loglik = sum(data$cv_loglik),
    converged = converged,
    pdHess = pdHess,
    max_gradients = max_grad
  )
}

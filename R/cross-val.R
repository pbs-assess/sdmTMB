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

dstudent <- function(x, df, mean, sd, ncp, log = FALSE) {
  # from metRology::dt.scaled()
  if (!log) {
    return(stats::dt((x - mean) / sd, df, ncp = ncp, log = FALSE) / sd)
  } else {
    return(stats::dt((x - mean) / sd, df, ncp = ncp, log = TRUE) - log(sd))
  }
}

ll_student <- function(object, withheld_y, withheld_mu) {
  .sd <- exp(object$model$par[["ln_phi"]])
  dstudent(x = withheld_y, df = object$tmb_data$df, mean = withheld_mu, sd = .sd, log = TRUE)
}

ll_nbinom1 <- function(object, withheld_y, withheld_mu) {
  phi <- exp(object$model$par[["ln_phi"]])
  stats::dnbinom(x = withheld_y, size = withheld_mu / phi, mu = withheld_mu, log = TRUE)
}

ll_nbinom2 <- function(object, withheld_y, withheld_mu) {
  phi <- exp(object$model$par[["ln_phi"]])
  stats::dnbinom(x = withheld_y, size = phi, mu = withheld_mu, log = TRUE)
}

# no longer used within sdmTMB_cv(); uses TMB report() instead
ll_sdmTMB <- function(object, withheld_y, withheld_mu) {
  family_func <- switch(object$family$family,
    gaussian = ll_gaussian,
    tweedie = ll_tweedie,
    binomial = ll_binomial,
    lognormal = ll_lognormal,
    student = ll_student,
    Gamma = ll_gamma,
    nbinom1 = ll_nbinom1,
    nbinom2 = ll_nbinom2,
    cli_abort(paste0(
      object$family$family, " not yet implemented. ",
      "Please file an issue on GitHub."
    ))
  )
  family_func(object, withheld_y, withheld_mu)
}

#' Cross validation with sdmTMB models
#'
#' Facilitates cross validation with sdmTMB models. Returns the log likelihood
#' of left-out data, which is similar in spirit to the ELPD (expected log
#' pointwise predictive density). The function has an option for
#' leave-future-out cross validation. By default, the function creates folds
#' randomly but folds can be manually assigned via the `fold_ids` argument.
#'
#' @param formula Model formula.
#' @param data A data frame.
#' @param mesh Output from [make_mesh()]. If supplied, the mesh will be constant
#'   across folds.
#' @param mesh_args Arguments for [make_mesh()]. If supplied, the mesh will be
#'   reconstructed for each fold.
#' @param time The name of the time column. Leave as `NULL` if this is only
#'   spatial data.
#' @param k_folds Number of folds.
#' @param fold_ids Optional vector containing user fold IDs. Can also be a
#'   single string, e.g. `"fold_id"` representing the name of the variable in
#'   `data`. Ignored if `lfo` is TRUE
#' @param lfo Whether to implement leave-future-out (LFO) cross validation where
#'   data are used to predict future folds. `time` argument in [sdmTMB()] must
#'   be specified. See Details section below.
#' @param lfo_forecast If `lfo = TRUE`, number of time steps to forecast. Time
#'   steps 1, ..., T are used to predict T + `lfo_forecast` and the last
#'   forecasted time step is used for validation. See Details section below.
#' @param lfo_validations If `lfo = TRUE`, number of times to step through the
#'   LFOCV process. Defaults to 5. See Details section below.
#' @param parallel If `TRUE` and a [future::plan()] is supplied, will be run in
#'   parallel.
#' @param use_initial_fit Fit the first fold and use those parameter values
#'   as starting values for subsequent folds? Can be faster with many folds.
#' @param future_globals A character vector of global variables used within
#'   arguments if an error is returned that \pkg{future.apply} can't find an
#'   object. This vector is appended to `TRUE` and passed to the argument
#'   `future.globals` in [future.apply::future_lapply()]. Useful if global
#'   objects are used to specify arguments like priors, families, etc.
#' @param spde **Depreciated.** Use `mesh` instead.
#' @param ... All other arguments required to run [sdmTMB()] model with the
#'   exception of `weights`, which are used to define the folds.
#'
#' @export
#' @return
#' A list:
#' * `data`: Original data plus columns for fold ID, CV predicted value,
#'           and CV log likelihood.
#' * `models`: A list of models; one per fold.
#' * `fold_loglik`: Sum of left-out log likelihoods per fold.
#' * `sum_loglik`: Sum of `fold_loglik` across all left-out data.
#' * `pdHess`: Logical vector: Hessian was invertible each fold?
#' * `converged`: Logical: all `pdHess` `TRUE`?
#' * `max_gradients`: Max gradient per fold.
#'
#' Prior to \pkg{sdmTMB} version '0.3.0.9002', `elpd` was incorrectly returned as
#' the log average likelihood, which is another metric you could compare models
#' with, but not ELPD. For maximum likelihood, [ELPD is equivalent in spirit to
#' the sum of the log likelihoods](https://github.com/pbs-assess/sdmTMB/issues/235).
#'
#' @details
#' **Parallel processing**
#'
#' Parallel processing can be used by setting a `future::plan()`.
#'
#' For example:
#'
#' ```
#' library(future)
#' plan(multisession)
#' # now use sdmTMB_cv() ...
#' ```
#'
#' **Leave-future-out cross validation (LFOCV)**
#'
#' An example of LFOCV with 9 time steps, `lfo_forecast = 1`, and
#' `lfo_validations = 2`:
#'
#' - Fit data to time steps 1 to 7, predict and validate step 8.
#' - Fit data to time steps 1 to 8, predict and validate step 9.
#'
#' An example of LFOCV with 9 time steps, `lfo_forecast = 2`, and
#' `lfo_validations = 3`:
#'
#' - Fit data to time steps 1 to 5, predict and validate step 7.
#' - Fit data to time steps 1 to 6, predict and validate step 8.
#' - Fit data to time steps 1 to 7, predict and validate step 9.
#'
#' See example below.
#'
#' @examples
#' mesh <- make_mesh(pcod, c("X", "Y"), cutoff = 25)
#'
#' # Set parallel processing first if desired with the future package.
#' # See the Details section above.
#'
#' m_cv <- sdmTMB_cv(
#'   density ~ 0 + depth_scaled + depth_scaled2,
#'   data = pcod, mesh = mesh,
#'   family = tweedie(link = "log"), k_folds = 2
#' )
#'
#' m_cv$fold_loglik
#' m_cv$sum_loglik
#'
#' head(m_cv$data)
#' m_cv$models[[1]]
#' m_cv$max_gradients
#'
#' \donttest{
#' # Create mesh each fold:
#' m_cv2 <- sdmTMB_cv(
#'   density ~ 0 + depth_scaled + depth_scaled2,
#'   data = pcod, mesh_args = list(xy_cols = c("X", "Y"), cutoff = 20),
#'   family = tweedie(link = "log"), k_folds = 2
#' )
#'
#' # Use fold_ids:
#' m_cv3 <- sdmTMB_cv(
#'   density ~ 0 + depth_scaled + depth_scaled2,
#'   data = pcod, mesh = mesh,
#'   family = tweedie(link = "log"),
#'   fold_ids = rep(seq(1, 3), nrow(pcod))[seq(1, nrow(pcod))]
#' )
#'
#' # LFOCV:
#' m_lfocv <- sdmTMB_cv(
#'   present ~ s(year, k = 4),
#'   data = pcod,
#'   mesh = mesh,
#'   lfo = TRUE,
#'   lfo_forecast = 2,
#'   lfo_validations = 3,
#'   family = binomial(),
#'   spatiotemporal = "off",
#'   time = "year" # must be specified
#' )
#'
#' # See how the LFOCV folds were assigned:
#' example_data <- m_lfocv$models[[1]]$data
#' table(example_data$cv_fold, example_data$year)
#' }
sdmTMB_cv <- function(
    formula, data, mesh_args, mesh = NULL, time = NULL,
    k_folds = 8, fold_ids = NULL,
    lfo = FALSE,
    lfo_forecast = 1,
    lfo_validations = 5,
    parallel = TRUE,
    use_initial_fit = FALSE,
    future_globals = NULL,
    spde = deprecated(),
    ...) {
  if (k_folds < 1) cli_abort("`k_folds` must be >= 1.")

  if (is_present(spde)) {
    deprecate_warn("0.0.20", "sdmTMB_cv(spde)", "sdmTMB_cv(mesh)")
  } else {
    spde <- mesh
  }
  data[["_sdm_order_"]] <- seq_len(nrow(data))
  constant_mesh <- missing(mesh_args)
  if (missing(mesh_args)) mesh_args <- NULL
  if (missing(spde)) spde <- NULL
  if (lfo) fold_ids <- NULL
  # add column of fold_ids stratified across time steps
  if (is.null(time)) {
    time <- "_sdmTMB_time"
    data[[time]] <- 0L
  }
  if (is.null(fold_ids)) {
    if (lfo) {
      if (length(unique(data[[time]])) < (lfo_validations + lfo_forecast)) {
        cli_abort("Not enough time steps for the desired validation period. Either decrease `lfo_validations` or add more data.")
      }
      # Create lfo_validations + 1 folds, ordered sequentially
      data$cv_fold <- 1
      t_validate <- sort(unique(data[[time]]), decreasing = TRUE)
      for (t in seq(1, lfo_validations)) {
        # fold id increasing order + forecast
        data$cv_fold[data[[time]] == t_validate[t]] <- lfo_validations - t + 1 + lfo_forecast
      }
      data$cv_fold <- as.numeric(as.factor(data$cv_fold))
    } else {
      dd <- lapply(split(data, data[[time]]), function(x) {
        x$cv_fold <- sample(rep(seq(1L, k_folds), nrow(x)), size = nrow(x))
        x
      })
      data <- do.call(rbind, dd)
    }
    fold_ids <- "cv_fold"
  } else {
    # fold_ids passed in; can be numeric, or a named column in `data`
    data$cv_fold <- NA
    if (length(fold_ids) == nrow(data)) {
      data$cv_fold <- fold_ids
    }
    if (length(fold_ids) == 1L && is.character(fold_ids)) {
      if (!fold_ids %in% names(data)) {
        cli_abort("Name of fold identifier not found in data.")
      }
      data$cv_fold <- data[[fold_ids]]
    }
    if (length(fold_ids) > 1 && length(fold_ids) < nrow(data)) {
      cli_abort("Dimension of `fold_ids` doesn't match data and is not a named variable.")
    }
    if (length(which(is.na(data$cv_fold))) > 0) {
      cli_abort("NAs found in `fold_ids`; please check `fold_ids` are specified correctly.")
    }
    k_folds <- length(unique(data$cv_fold))
  }
  if (time == "_sdmTMB_time") { # undo changes above, make time NULL
    data[["_sdmTMB_time"]] <- NULL
    time <- NULL
  }

  dot_args <- as.list(substitute(list(...)))[-1L]
  if ("weights" %in% names(dot_args)) {
    cli_abort("`weights` cannot be specified within sdmTMB_cv().")
  }
  if ("offset" %in% names(dot_args)) {
    .offset <- eval(dot_args$offset)
    if (parallel && !is.character(.offset) && !is.null(.offset)) {
      cli_abort("We recommend using a character value for 'offset' (indicating the column name) when applying parallel cross validation.")
    }
  } else {
    .offset <- NULL
  }

  if (k_folds > 1) {
    # data in kth fold get weight of 0:
    weights <- ifelse(data$cv_fold == 1L, 0, 1)
  } else {
    weights <- rep(1, nrow(data))
  }
  if (lfo) weights <- ifelse(data$cv_fold == 1L, 1, 0)

  if (use_initial_fit) {
    # run model on first fold to get starting values:

    if (!constant_mesh) {
      if (lfo) {
        dat_fit <- data[data$cv_fold == 1L, , drop = FALSE]
      } else {
        dat_fit <- data[data$cv_fold != 1L, , drop = FALSE]
      }

      mesh_args[["data"]] <- dat_fit
      mesh <- do.call(make_mesh, mesh_args)
    } else {
      mesh <- spde
      dat_fit <- data
    }
    dot_args <- list(dot_args)[[1]]
    dot_args$offset <- NULL
    .args <- c(list(
      data = dat_fit, formula = formula, time = time, mesh = mesh,
      weights = weights, offset = .offset
    ), dot_args)
    fit1 <- do.call(sdmTMB, .args)
  }

  fit_func <- function(k) {
    if (lfo) {
      weights <- ifelse(data$cv_fold <= k, 1, 0)
    } else {
      # data in kth fold get weight of 0:
      weights <- ifelse(data$cv_fold == k, 0, 1)
    }

    if (k == 1L && use_initial_fit) {
      object <- fit1
    } else {
      if (!constant_mesh) {
        if (lfo) {
          dat_fit <- data[data$cv_fold <= k, , drop = FALSE]
        } else {
          dat_fit <- data[data$cv_fold != k, , drop = FALSE]
        }
        mesh_args[["data"]] <- dat_fit
        mesh <- do.call(make_mesh, mesh_args)
      } else {
        mesh <- spde
        dat_fit <- data
      }
      dot_args <- as.list(substitute(list(...)))[-1L] # re-evaluate here! issue #54
      dot_args <- list(...)
      dot_args$offset <- NULL
      args <- c(list(
        data = dat_fit, formula = formula, time = time, mesh = mesh, offset = .offset,
        weights = weights, previous_fit = if (use_initial_fit) fit1 else NULL
      ), dot_args)
      object <- do.call(sdmTMB, args)
      # if (max(object$gradients) > 0.01) {
      # object <- run_extra_optimization(object, nlminb_loops = 1L, newton_loops = 0L)
      # }
    }

    if (lfo) {
      cv_data <- data[data$cv_fold == (k + lfo_forecast), , drop = FALSE]
    } else {
      cv_data <- data[data$cv_fold == k, , drop = FALSE]
    }

    # FIXME: only use TMB report() below to be faster!
    # predict for withheld data:
    predicted <- predict(object, newdata = cv_data, type = "response")
    cv_data$cv_predicted <- predicted$est
    response <- get_response(object$formula[[1]])
    withheld_y <- predicted[[response]]
    withheld_mu <- cv_data$cv_predicted

    # FIXME: get LFO working with the TMB report() option below!
    # calculate log likelihood for each withheld observation:
    # trickery to get the log likelihood of the withheld data directly
    # from the TMB report():
    if (!lfo) {
      tmb_data <- object$tmb_data
      tmb_data$weights_i <- ifelse(tmb_data$weights_i == 1, 0, 1) # reversed
      new_tmb_obj <- TMB::MakeADFun(
        data = tmb_data,
        parameters = get_pars(object),
        map = object$tmb_map,
        random = object$tmb_random,
        DLL = "sdmTMB",
        silent = TRUE
      )
      lp <- object$tmb_obj$env$last.par.best
      r <- new_tmb_obj$report(lp)
      cv_loglik <- -1 * r$jnll_obs
      cv_data$cv_loglik <- cv_loglik[tmb_data$weights_i == 1]
    } else { # old method; doesn't work with delta models!
      cv_data$cv_loglik <- ll_sdmTMB(object, withheld_y, withheld_mu)
    }

    ## test
    # x2 <- ll_sdmTMB(object, withheld_y, withheld_mu)
    # identical(round(cv_data$cv_loglik, 6), round(x2, 6))
    # cv_data$cv_loglik <- ll_sdmTMB(object, withheld_y, withheld_mu)

    list(
      data = cv_data,
      model = object,
      pdHess = object$sd_report$pdHess,
      max_gradient = max(abs(object$gradients)),
      bad_eig = object$bad_eig
    )
  }

  if (requireNamespace("future.apply", quietly = TRUE) && parallel) {
    message(
      "Running fits with `future.apply()`.\n",
      "Set a parallel `future::plan()` to use parallel processing."
    )
    if (!is.null(future_globals)) {
      fg <- structure(TRUE, add = future_globals)
    } else {
      fg <- TRUE
    }
    if (lfo) {
      out <- future.apply::future_lapply(seq_len(lfo_validations), fit_func, future.seed = TRUE, future.globals = fg)
    } else {
      out <- future.apply::future_lapply(seq_len(k_folds), fit_func, future.seed = TRUE, future.globals = fg)
    }
  } else {
    message(
      "Running fits sequentially.\n",
      "Install the future and future.apply packages,\n",
      "set a parallel `future::plan()`, and set `parallel = TRUE` to use parallel processing."
    )
    if (lfo) {
      out <- lapply(seq_len(lfo_validations), fit_func)
    } else {
      out <- lapply(seq_len(k_folds), fit_func)
    }
  }

  models <- lapply(out, `[[`, "model")
  data <- lapply(out, `[[`, "data")
  fold_cv_ll <- vapply(data, function(.x) sum(.x$cv_loglik), FUN.VALUE = numeric(1L))
  data <- do.call(rbind, data)
  data <- data[order(data[["_sdm_order_"]]), , drop = FALSE]
  data[["_sdm_order_"]] <- NULL
  data[["_sdmTMB_time"]] <- NULL
  row.names(data) <- NULL
  pdHess <- vapply(out, `[[`, "pdHess", FUN.VALUE = logical(1L))
  max_grad <- vapply(out, `[[`, "max_gradient", FUN.VALUE = numeric(1L))
  converged <- all(pdHess)
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

log_sum_exp <- function(x) {
  max_x <- max(x)
  max_x + log(sum(exp(x - max_x)))
}

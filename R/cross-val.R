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

ll_sdmTMB <- function(object, withheld_y, withheld_mu) {
  family_func <- switch(object$family$family,
    gaussian = ll_gaussian,
    tweedie = ll_tweedie,
    binomial = ll_binomial
  )
  family_func(object, withheld_y, withheld_mu)
}

#' Save log likelihoods of k-fold cross-validation for sdmTMB models
#'
#' @param formula Model formula.
#' @param data A data frame.
#' @param time The name of the time column. Leave as `NULL` if this is only spatial data.
#' @param x Name of the column with X coordinates.
#' @param y Name of the column with Y coordinates.
#' @param k_folds Number of folds.
#' @param fold_ids Optional input name of column containing user chosen fold
#'   ids.
#' @param n_knots The number of knots.
#' @param spde_function A function that takes 3 arguments (`x`, `y`, `n_knots`)
#'   (and `mesh` if `knot_type = "fixed"`) and returns a list structure that
#'   matches output of [make_spde()].
#' @param knot_type Should the mesh knots be fixed across each fold (`"fixed"`)
#'   or generated separately for each fold (`"unique"`)?
#' @param seed A seed to ensure the fold pattern is the same for comparison
#'   purposes.
#' @param ... All other arguments required to run sdmTMB model with the
#'   exception of `data` and `spde` which are redefined for each fold within the
#'   function.
#'
#' @export
#'
#' @examples
#' d <- subset(pcod, year >= 2011) # subset for example speed
#' x <- sdmTMB_cv(
#'   formula = density ~ 0 + depth_scaled + depth_scaled2 + as.factor(year),
#'   d, family = tweedie(link = "log"), time = "year", x = "X", y = "Y",
#'   n_knots = 30, k_folds = 2
#' )
#' x$sum_loglik
#' str(x$data)
#' x$models[[1]]
#' x$models[[2]]
#'
#' # proof the knots are the same but the data change if `knot_type = "fixed"`:
#' plot_spde(x$models[[1]]$spde)
#' plot_spde(x$models[[2]]$spde)

sdmTMB_cv <- function(formula, data, x, y, time = NULL,
                      k_folds = 10, fold_ids = NULL, n_knots = NULL,
                      spde_function = make_spde,
                      knot_type = c("fixed", "unique"),
                      seed = 999, ...) {
  set.seed(seed)
  data[["_sdm_order_"]] <- seq_len(nrow(data))
  knot_type <- match.arg(knot_type)

  if (is.null(time)) {
    time <- "_sdmTMB_time"
    data[[time]] <- 0L
  }

  # add column of fold_ids stratified across time steps
  if (is.null(fold_ids)) {
    dd <- lapply(split(data, data[[time]]), function(x) {
      x$cv_fold <- sample(rep(seq(1L, k_folds), nrow(x)), size = nrow(x))
      x
    })
    data <- do.call(rbind, dd)
    fold_ids <- "cv_fold"
  } else {
    data$cv_fold <- fold_ids
  }

  if (identical(knot_type, "fixed")) {
    spde_global <- spde_function(data[[x]], data[[y]], n_knots = n_knots)
  }

  out <- lapply(seq_len(k_folds), function(k) {
    if (k_folds > 1) {
      d_fit <- data[data[[fold_ids]] != k, , drop = FALSE]
      d_withheld <- data[data[[fold_ids]] == k, , drop = FALSE]
    } else {
      d_fit <- data
      d_withheld <- data
    }

    # build mesh for training data
    if (identical(knot_type, "fixed")) {
      d_fit_spde <- spde_function(d_fit[[x]], d_fit[[y]], n_knots = n_knots,
        mesh = spde_global$mesh)
    } else {
      d_fit_spde <- spde_function(d_fit[[x]], d_fit[[y]], n_knots = n_knots)
    }

    # run model
    object <- sdmTMB(data = d_fit, formula = formula, time = time, spde = d_fit_spde, ...)

    # predict for withheld data
    predicted <- predict(object, newdata = d_withheld, xy_cols = c(x, y))
    cv_data <- d_withheld
    cv_data$cv_predicted <- object$family$linkinv(predicted$est)
    response <- get_response(object$formula)
    withheld_y <- predicted[[response]]
    withheld_mu <- cv_data$cv_predicted

    # calculate log likelihood for each withheld observationn
    cv_data$cv_loglik <- ll_sdmTMB(object, withheld_y, withheld_mu)

    list(
      data = cv_data,
      model = object,
      pdHess = object$sd_report$pdHess,
      max_gradient = max(abs(object$gradients)),
      bad_eig = object$bad_eig
    )
  })
  models <- lapply(out, `[[`, "model")
  data <- lapply(out, `[[`, "data")
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
    sum_loglik = sum(data$cv_loglik),
    converged = converged,
    max_grad = max(max_grad)
  )
}

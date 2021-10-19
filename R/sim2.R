#' Simulate from a spatial/spatiotemporal model
#'
#' This version (vs. [sdmTMB_sim()]) is set up to take a formula and a data
#' frame and is easier to use if you want different spatial observations (and
#' covariates) for each time slice. Eventually [sdmTMB_sim()] will be
#' depreciated in favour of this version.
#'
#' @param formula A *one-sided* formula describing the fixed-effect structure.
#' @param data A data frame containing the predictors described in `formula` and the
#'   time column if `time` is specified.
#' @param spde Output from [make_mesh()].
#' @param time The time column name.
#' @param family Family as in [sdmTMB()].
#' @param range Parameter that controls the decay of spatial correlation.
#' @param X An optional covariate design matrix. If omitted and `B` is not
#'   `NULL`, will be set to standard normal draws.
#' @param B A vector of beta values (fixed-effect coefficient values).
#' @param rho Spatiotemporal correlation between years; should be between -1 and
#'   1.
#' @param sigma_O SD of spatial process (Omega).
#' @param sigma_E SD of spatiotemporal process (Epsilon). Can be scalar or
#'   vector for time-varying model.
#' @param phi Observation error scale parameter (e.g., SD in Gaussian).
#' @param tweedie_p Tweedie p (power) parameter; between 1 and 2.
#' @param df Student-t degrees of freedom.
#' @param seed A value with which to set the random seed.
#' @param size Specific for the binomial family, vector representing binomial N.
#'   If not included, defaults to 1 (Bernoulli)
#'
#' @return A data frame where:
#' * The 1st column is the time variable (if present).
#' * The 2nd and 3rd columns are the spatial coordinates.
#' * `omega_s` represents the simulated spatial random effects.
#' * `epsilon_st` represents the simulated spatiotemporal random effects.
#' * `eta` is the true value in link space
#' * `mu` is the true value in inverse link space.
#' * `observed` represents the simulated process with observation error.
#' * The remaining columns are the fixed-effect model matrix.
#' @export
#'
#' @examples
#' if (inla_installed()) {
#'
#'   set.seed(123)
#'   # a1 is a fake predictor:
#'   predictor_dat <- data.frame(
#'     X = runif(300), Y = runif(300),
#'     a1 = rnorm(300), year = rep(1:6, each = 50)
#'   )
#'   mesh <- make_mesh(predictor_dat, xy_cols = c("X", "Y"), cutoff = 0.1)
#'
#'   sim_dat <- sdmTMB_sim2(
#'     formula = ~ 1 + a1,
#'     data = predictor_dat,
#'     time = "year",
#'     spde = mesh,
#'     family = gaussian(link = "identity"),
#'     range = 0.5,
#'     sigma_E = 0.1,
#'     phi = 0.1,
#'     sigma_O = 0.2,
#'     seed = 3542,
#'     B = c(0.2, -0.4) # B0 = intercept, B1 = a1 slope
#'   )
#'
#'   fit <- sdmTMB(observed ~ a1, data = sim_dat, spde = mesh, time = "year")
#'   fit
#'
#'   if (require("ggplot2", quietly = TRUE)) {
#'     ggplot(sim_dat, aes(X, Y, colour = observed)) +
#'       geom_point() +
#'       facet_wrap(~year) +
#'       scale_color_gradient2()
#'   }
#' }
sdmTMB_sim2 <- function(formula,
                        data,
                        spde,
                        time = NULL,
                        family = gaussian(link = "identity"),
                        # time_varying = NULL,
                        range,
                        B = 0,
                        rho = 0,
                        sigma_O = 0.1,
                        sigma_E = 0,
                        phi = 0.1,
                        tweedie_p = 1.5,
                        df = 3,
                        size = NULL,
                        seed = sample.int(1e6, 1)) {
  mesh <- spde
  betas <- B
  if (!requireNamespace("INLA", quietly = TRUE)) {
    stop("INLA must be installed to use this function.", call. = FALSE)
  }
  assert_that(class(mesh) %in% c("inla.mesh", "sdmTMBmesh"))
  assert_that(tweedie_p > 1, tweedie_p < 2)
  assert_that(df >= 1)
  assert_that(range > 0)
  assert_that(rho >= -1, rho <= 1)
  assert_that(sigma_O >= 0, all(sigma_E >= 0), phi > 0)

  if (class(mesh) == "sdmTMBmesh") {
    mesh <- mesh$mesh
  }

  response <- get_response(formula)
  if (length(response) == 0L) {
    formula <- as.formula(paste("sdmTMB_response_", paste(as.character(formula), collapse = "")))
    data[["sdmTMB_response_"]] <- 0.1 # fake! does nothing but lets sdmTMB parse the formula
  }
  # get tmb_data structure; parsed model matrices etc.:
  fit <- sdmTMB(
    formula = formula, data = data, spde = spde, time = time,
    family = family, time_varying = NULL, do_fit = FALSE
  )

  X_ij <- fit$tmb_data$X_ij
  n_t <- fit$tmb_data$n_t

  n_covariates <- length(B)
  assert_that(ncol(X_ij) == length(B),
    msg = "Number of B parameters does not match model matrix columns implied by the formula."
  )

  coords <- spde$loc_xy

  rspde_attr_O <- get_rspde3_attributes(
    coords = coords, sigma = sigma_O, range = range,
    mesh = mesh
  )

  if (sigma_O > 0) {
    omega_s <- rspde3(rspde_attr_O, seed = seed)
  } else {
    omega_s <- rep(0, length(x))
  }

  # test whether sigma_E_zero
  if (length(sigma_E) %in% c(1L, n_t) == FALSE) {
    stop("Error: sigma_E must be a scalar or of length time_steps", call. = FALSE)
  }
  if (length(sigma_E) == 1L) {
    rspde_attr_E <- get_rspde3_attributes(
      coords = coords, sigma = sigma_E, range = range,
      mesh = mesh
    )
    sigma_E <- rep(sigma_E, n_t)
    rspde_attr_E <- rep(list(rspde_attr_E), n_t)
  } else {
    rspde_attr_E <- lapply(sigma_E, function(.x) {
      get_rspde3_attributes(
        coords = coords,
        sigma = .x, range = range,
        mesh = mesh
      )
    })
  }

  epsilon_st <- list() # spatiotemporal random effects
  if (sigma_E[[1]] > 0 && n_t > 1L) {
    for (i in seq_len(n_t)) {
      if (i == 1 || rho == 0) {
        epsilon_st[[i]] <- rspde3(rspde_attr_E[[i]], seed = seed * i)
      } else { # AR1 and not first time slice:
        epsilon_st[[i]] <- rho * epsilon_st[[i - 1]] +
          sqrt(1 - rho^2) * rspde3(rspde_attr_E[[i]], seed = seed * i)
      }
    }
  } else {
    epsilon_st <- list(rep(0, nrow(X_ij)))
  }

  if (!is.null(time)) {
    # only retain observations each time slice:
    time_vec <- sort(unique(data[[time]]))
    for (i in seq_len(n_t)) {
      epsilon_st[[i]] <- epsilon_st[[i]][which(data[[time]] == time_vec[[1]])]
    }
  }
  epsilon_st <- do.call("c", epsilon_st)

  eta <- X_ij %*% B

  d <- list()
  if (!is.null(time)) d[[time]] <- data[[time]]
  d[[spde$xy_cols[1]]] <- coords[, 1]
  d[[spde$xy_cols[2]]] <- coords[, 2]
  d[["omega_s"]] <- omega_s
  d[["epsilon_st"]] <- epsilon_st
  d[["eta"]] <- eta
  d <- do.call("data.frame", d)
  d$eta <- d$eta + d$omega_s + d$epsilon_st

  d$mu <- do.call(family$linkinv, list(d$eta))
  N <- nrow(d)

  if (is.null(size)) size <- rep(1, N)
  d$observed <- switch(family$family,
    gaussian  = stats::rnorm(N, mean = d$mu, sd = phi),
    binomial  = stats::rbinom(N, size = size, prob = d$mu),
    tweedie   = fishMod::rTweedie(N, mu = d$mu, phi = phi, p = tweedie_p),
    Beta      = stats::rbeta(N, d$mu * phi, (1 - d$mu) * phi),
    Gamma     = stats::rgamma(N, shape = phi, scale = d$mu / phi),
    nbinom2   = stats::rnbinom(N, size = phi, mu = d$mu),
    poisson   = stats::rpois(N, lambda = d$mu),
    student   = rstudent(N, d$mu, sigma = phi, nu = df),
    lognormal = stats::rlnorm(N, meanlog = log(d$mu) - (phi^2) / 2, sdlog = phi),
    stop("Family not found.", call. = FALSE)
  )

  cbind(d, X_ij)
}

# modified from https://becarioprecario.bitbucket.io/spde-gitbook/ch-spacetime.html
get_rspde3_attributes <- function(coords, mesh, sigma = 1, range, variance = sigma^2, alpha = 2,
                                  kappa = sqrt(8 * (alpha - 1)) / range) {
  theta <- c(-0.5 * log(4 * pi * variance * kappa^2), log(kappa))
  if (missing(mesh)) {
    stop("`mesh` must be specified.", call. = FALSE)
  } else {
    attributes <- list(mesh = mesh)
  }
  attributes$spde <- INLA::inla.spde2.matern(attributes$mesh, alpha = alpha)
  attributes$Q <- INLA::inla.spde2.precision(attributes$spde, theta = theta)
  attributes$A <- INLA::inla.mesh.project(mesh = attributes$mesh, loc = coords)$A
  attributes
}

# modified from https://becarioprecario.bitbucket.io/spde-gitbook/ch-spacetime.html
rspde3 <- function(attributes = NULL, seed = 0L, n = 1L) {
  if (as.integer(n) == 1L) {
    result <- drop(
      attributes$A %*% INLA::inla.qsample(
        Q = attributes$Q,
        seed = if (missing(seed)) 0L else seed,
        num.threads = 1L,
        constr = attributes$spde$f$extraconstr
      )
    )
  } else {
    stop("`n` must be 1.", call. = FALSE)
  }
  result <- drop(result)
  result <- as.matrix(result)
  colnames(result) <- NULL
  result
}

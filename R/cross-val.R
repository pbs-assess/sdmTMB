#' Simulate from a spatial/spatiotemporal model
#'
#' @param mesh Output from [make_mesh()] or a mesh directly from INLA.
#' @param x A vector of x coordinates. Should match `mesh`.
#' @param y A vector of y coordinates. Should match `mesh`.
#' @param range Parameter that controls the decay of spatial correlation.
#' @param X An optional covariate design matrix formatted as a list with each
#'   element of the list representing a slice in time. If ommitted and `betas`
#'   is not `NULL`, will be set to standard normal.
#' @param betas A vector of beta values (design-matrix fixed-effect coefficient
#'   values). If a random walk (`sigma_V > 0`), these are the starting values.
#' @param family Family as in [sdmTMB()].
#' @param time_steps The number of time steps.
#' @param rho Spatiotemporal correlation between years;
#'   should be between -1 and 1.
#' @param sigma_O SD of spatial process (Omega).
#' @param sigma_E SD of spatiotemporal process (Epsilon).
#' @param sigma_V A vector of standard deviations of time-varying random walk on
#'   parameters. Set to 0 for parameters that should not vary through time.
#' @param phi Observation error scale parameter.
#' @param thetaf Tweedie p (power) parameter; between 1 and 2.
#' @param df Student-t degrees of freedom.
#' @param seed A value with which to set the random seed.
#' @param list Logical for whether output is in list format. If `TRUE`,
#'    data is in list element 1 and input values in element 2.
#'
#' @return A data frame where:
#' * `omega_s` represents the simulated spatial random effects.
#' * `epsilon_st` represents the simulated spatiotemporal random effects.
#' * `eta` is the true value in link space
#' * `mu` is the true value in inverse link space.
#' * `observed` represents the simulated process with observation error.
#' * `b_...` contain the beta values for each covariate used to simulate each time slice.
#' * `cov_...` covariate values for each observation.

#' @export
#'
#' @examples
#' \donttest{
#' set.seed(42)
#' x <- runif(50, -1, 1)
#' y <- runif(50, -1, 1)
#' N <- length(x)
#' time_steps <- 6
#' X <- model.matrix(~ x1, data.frame(x1 = rnorm(N * time_steps)))
#' loc <- data.frame(x = x, y = y)
#' mesh <- make_mesh(loc, xy_cols = c("x", "y"), cutoff = 0.1)
#' s <- sdmTMB_sim(
#'   x = x, y = y, mesh = mesh, X = X,
#'   betas = c(0.5, 0.7), time_steps = time_steps, rho = 0.5,
#'   phi = 0.2, range = 0.8, sigma_O = 0, sigma_E = 0.3,
#'   seed = 123, family = gaussian()
#' )
#'
#' mesh <- make_mesh(s, xy_cols = c("x", "y"), cutoff = 0.1)
#' m <- sdmTMB(
#'   data = s, formula = observed ~ x1,
#'   time = "time", spde = mesh,
#'   ar1_fields = TRUE, include_spatial = FALSE
#' )
#' tidy(m, conf.int = TRUE)
#' tidy(m, "ran_pars", conf.int = TRUE)
#' }

sdmTMB_sim <- function(mesh,
  x,
  y,
  range,
  time_steps = 1L,
  X = NULL,
  betas = NULL,
  family = gaussian(link = "identity"),
  rho = 0,
  sigma_O = 0.1,
  sigma_E = 0,
  sigma_V = rep(0, length(betas)),
  phi = 0.01,
  thetaf = 1.5,
  df = 3,
  seed = sample.int(1e6, 1),
  list = FALSE) {
  assert_that(is.numeric(x), is.numeric(y))
  assert_that(is.null(dim(x)), is.null(dim(y)))
  assert_that(identical(length(x), length((y))))
  assert_that(class(mesh) %in% c("inla.mesh", "sdmTMBmesh"))
  assert_that(thetaf > 1, thetaf < 2)
  assert_that(df >= 1)
  assert_that(time_steps >= 1)
  assert_that(range > 0)
  assert_that(rho >= -1, rho <= 1)
  assert_that(sigma_O >= 0, sigma_E >= 0, all(sigma_V >= 0), phi > 0)
  if (!is.null(X)) assert_that(!is.null(betas))
  if (!is.null(betas) && !is.null(X)) assert_that(ncol(X) == length(betas))
  assert_that(length(betas) == length(sigma_V))
  if (!is.null(X)) assert_that(time_steps * length(x) == nrow(X))

  if (class(mesh) == "sdmTMBmesh") {
    mesh <- mesh$mesh
  }
  n_covariates <- length(betas)
  coords <- cbind(x, y)
  if (sigma_O > 0) {
    omega_s <- rspde2(coords, sigma = sigma_O, range = range, mesh = mesh, seed = seed)
  } else {
    omega_s <- rep(0, length(x))
  }
  epsilon_st <- list() # spatiotemporal random effects

  if (time_steps > 1L && sigma_E > 0) {
    for (i in seq_len(time_steps)) {
      if (i == 1 || rho == 0) {
        epsilon_st[[i]] <-
          rspde2(coords, sigma = sigma_E, range = range, mesh = mesh, seed = seed * i)
      } else { # AR1 and not first time slice:
        epsilon_st[[i]] <- rho * epsilon_st[[i - 1]] +
          sqrt(1 - rho^2) * # stationary AR1
          rspde2(coords, sigma = sigma_E, range = range, mesh = mesh, seed = seed * i)
      }
    }
  } else {
    epsilon_st <- list(rep(0, length(x)))
  }
  epsilon_st <- do.call("c", epsilon_st)

  # create betas for each covariate k
  if (n_covariates > 0) {
    b <- matrix(ncol = n_covariates, nrow = time_steps)
    b[1, ] <- betas
    if (time_steps > 1) {
      for (k in seq_len(n_covariates)) {
        for (i in seq(2, time_steps)) {
          b[i, k] <- b[i - 1, k] + stats::rnorm(1, 0, sigma_V[k])
        }
      }
    }
    if (is.null(X)) {
      cov_mat <- matrix(stats::rnorm(length(x) * n_covariates * time_steps, sd = 1),
        ncol = n_covariates, nrow = length(x) * time_steps
      )
      cov_mat <- as.data.frame(cov_mat)
      names(cov_mat) <- gsub("V", "cov", names(cov_mat))
    } else {
      cov_mat <- X
    }

    V <- b[rep(seq_len(time_steps), each = length(x)), ]
    B <- as.data.frame(V)
    names(B) <- gsub("V", "b", names(B))
    cov_values <- B * cov_mat
    eta <- rowSums(cov_values)
  } else {
    # empty vector for eta in simulations with no fixed effects
    eta <- vector("numeric", length = length(x) * time_steps)
  }

  d <- data.frame(
    time = rep(seq_len(time_steps), each = length(x)),
    x, y,
    omega_s = rep(omega_s, time_steps), epsilon_st = epsilon_st,
    eta = eta
  )
  d$eta <- d$eta + d$omega_s + d$epsilon_st

  d$mu <- do.call(family$linkinv, list(d$eta))
  N <- nrow(d)

  d$observed <- switch(family$family,
    gaussian  = stats::rnorm(N, mean = d$mu, sd = phi),
    binomial  = stats::rbinom(N, size = 1L, prob = d$mu),
    tweedie   = fishMod::rTweedie(N, mu = d$mu, phi = phi, p = thetaf),
    Beta      = stats::rbeta(N, d$mu * phi, (1 - d$mu) * phi),
    Gamma     = stats::rgamma(N, shape = 1 / (phi^2), scale = d$mu / (1 / (phi^2))),
    nbinom2   = stats::rnbinom(N, size = phi, mu = d$mu),
    poisson   = stats::rpois(N, lambda = d$mu),
    student   = rstudent(N, d$mu, sigma = phi, nu = df),
    lognormal = stats::rlnorm(N, meanlog = log(d$mu) - (phi^2) / 2, sdlog = phi),
    stop("Family not found.", call. = FALSE)
  )

  if (n_covariates > 0) {
    d <- cbind(d, B)
    d <- cbind(d, cov_mat)
  }
  if (list) {
    sim_out <- list(data = d, inputs = list(
      x = x,
      y = y,
      X = X,
      betas = betas,
      sigma_O = sigma_O,
      sigma_E = sigma_E,
      sigma_V = sigma_V,
      range = range,
      rho = rho,
      phi = phi,
      df = df,
      thetaf = thetaf
    ))
    sim_out
  } else {
    d
  }
}

# modified from https://becarioprecario.bitbucket.io/spde-gitbook/ch-spacetime.html
rspde2 <- function(coords, mesh, sigma = 1, range, variance = sigma^2, alpha = 2,
  kappa = sqrt(8 * (alpha - 1)) / range, n = 1,
  seed = 0L, return.attributes = FALSE) {
  theta <- c(-0.5 * log(4 * pi * variance * kappa^2), log(kappa))
  if (missing(mesh)) {
    stop("`mesh` must be specified.", call. = FALSE)
  }
  else {
    attributes <- list(mesh = mesh)
  }
  attributes$spde <- INLA::inla.spde2.matern(attributes$mesh, alpha = alpha)
  attributes$Q <- INLA::inla.spde2.precision(attributes$spde, theta = theta)
  attributes$A <- INLA::inla.mesh.project(mesh = attributes$mesh, loc = coords)$A
  # attributes$A <-  INLA::inla.spde.make.A(attributes$mesh, loc = coords)
  if (as.integer(n) == 1L) {
    result <- drop(
      attributes$A %*% INLA::inla.qsample(
        Q = attributes$Q,
        seed = if (missing(seed)) 0L else seed,
        # seed = 0L,
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

rstudent <- function(n, mu = 0, sigma = 1, nu = 3) {
  mu <- rep(mu, len = n)
  sigma <- rep(sigma, len = n)
  nu <- rep(nu, len = n)
  y <- stats::rnorm(n)
  z <- stats::rchisq(n, nu)
  mu + sigma * y * sqrt(nu / z)
}

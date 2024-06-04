## TODO: probably don't need both nproj and newdata; can drop nproj?

#' Project from an \pkg{sdmTMB} model using simulation
#'
#' @description `r lifecycle::badge("experimental")`
#'
#' @description The function enables projecting forward in time from an
#' \pkg{sdmTMB} model using a simulation approach for computational efficiency.
#' This can be helpful for calculating predictive intervals for long
#' projections where including those time elements in `extra_time` during model
#' estimation can be slow.
#'
#' @description Inspiration for this approach comes from the \pkg{VAST} function
#' `project_model()` by J.T. Thorson.
#'
#' @param object A fitted model from [sdmTMB()].
#' @param newdata A new data frame to predict on. Should contain all new time
#'   elements.
#' @param nproj Number of years to project.
#' @param nsim Number of simulations.
#' @param silent Silent?
#' @param sims_var Element to extract from the \pkg{TMB} report. Also see
#'   `return_tmb_report`.
#' @param model Linear predictor number to extract. Also see
#'   `return_tmb_report`.
#' @param sim_re A vector of `0`s and `1`s representing which random effects to
#'   simulate in the projection. Generally, leave this untouched. Order is:
#'   spatial fields, spatiotemporal fields, spatially varying coefficient
#'   fields, random intercepts, time-varying coefficients, smoothers.
#'   The default is to simulate spatiotemporal fields and time-varying
#'   coefficients, if present.
#' @param return_tmb_report Return the \pkg{TMB} report from `simulate()`? This
#'   lets you parse out whatever elements you want from the simulation including
#'   grabbing multiple elements from one set of simulations. See examples.
#' @param ... Passed to [predict.sdmTMB()].
#'
#' @references `project_model()` in the \pkg{VAST} package.
#'
#' @return
#' If `return_tmb_report = FALSE` (default): a matrix with N rows equal to the
#' number of rows in `newdata` and N columns equal to `nsim`.
#'
#' If `return_tmb_report = TRUE`: a list of TMB reports from `simulate()`. Run
#' `names()` on the output to see the options.
#'
#' TODO: fill in more details.
#' @export
#'
#' @examplesIf ggplot2_installed()
#' \donttest{
#' library(ggplot2)
#'
#' mesh <- make_mesh(dogfish, c("X", "Y"), cutoff = 25)
#' historical_years <- 2004:2022
#' to_project <- 10
#' future_years <- seq(max(historical_years) + 1, max(historical_years) + to_project)
#' all_years <- c(historical_years, future_years)
#' proj_grid <- replicate_df(wcvi_grid, "year", all_years)
#'
#' # we could fit our model like this, but for long projections, this becomes slow:
#' fit <- sdmTMB(
#'   catch_weight ~ 1,
#'   time = "year",
#'   offset = log(dogfish$area_swept),
#'   extra_time = all_years, #< note that all years here
#'   spatial = "on",
#'   spatiotemporal = "ar1",
#'   data = dogfish,
#'   mesh = mesh,
#'   family = tweedie(link = "log")
#' )
#'
#' # instead, we could fit our model like this and then take simulation draws
#' # from the projection time period:
#' fit2 <- sdmTMB(
#'   catch_weight ~ 1,
#'   time = "year",
#'   offset = log(dogfish$area_swept),
#'   extra_time = historical_years, #< does *not* include projection years
#'   spatial = "on",
#'   spatiotemporal = "ar1",
#'   data = dogfish,
#'   mesh = mesh,
#'   family = tweedie(link = "log")
#' )
#'
#' # we will only use 20 `nsim` so this example runs quickly
#' # you will likely want many more (> 200) in practice so the result
#' # is relatively stable
#'
#' set.seed(1)
#' out <- project(fit2, newdata = proj_grid, nproj = to_project, nsim = 20)
#' est_mean <- apply(out, 1, mean) # summarize however you'd like
#' est_se <- apply(out, 1, sd)
#'
#' # visualize:
#' proj_grid$est_mean <- est_mean
#' ggplot(subset(proj_grid, year > 2021), aes(X, Y, fill = est_mean)) +
#'   geom_raster() +
#'   facet_wrap(~year) +
#'   coord_fixed() +
#'   scale_fill_viridis_c() +
#'   ggtitle("Projection simulation (mean)")
#'
#' # if instead we wanted to grab, say, the spatiotemporal random field values,
#' # we can return the report and work with the raw output ourselves:
#'
#' set.seed(1)
#' out <- project(
#'   fit2,
#'   newdata = proj_grid, nproj = to_project,
#'   nsim = 20, # increase this
#'   return_tmb_report = TRUE #< difference from above example
#' )
#'
#' # here are the elements we could extract:
#' names(out[[1]])
#'
#' # 'epsilon_st_A_vec' are the 'epsilon_st' at every location in 'newdata':
#' eps <- lapply(out, \(x) x[["epsilon_st_A_vec"]][, 1])
#' eps <- do.call(cbind, eps)
#' eps_mean <- apply(eps, 1, mean) # summarize however you'd like
#' eps_se <- apply(eps, 1, sd)
#'
#' proj_grid$eps_mean <- eps_mean
#' ggplot(subset(proj_grid, year > 2021), aes(X, Y, fill = eps_mean)) +
#'   geom_raster() +
#'   facet_wrap(~year) +
#'   scale_fill_gradient2() +
#'   coord_fixed() +
#'   ggtitle("Projection simulation\n(spatiotemporal fields)")
#'
#' proj_grid$eps_se <- eps_se
#' ggplot(subset(proj_grid, year > 2021), aes(X, Y, fill = eps_se)) +
#'   geom_raster() +
#'   facet_wrap(~year) +
#'   scale_fill_viridis_c() +
#'   coord_fixed() +
#'   ggtitle("Projection simulation\n(spatiotemporal fields standard error)")
#' }

project <- function(
    object,
    newdata,
    nproj = 1,
    nsim = 1,
    silent = FALSE,
    sims_var = "eta_i",
    model = 1,
    sim_re = c(0, 1, 0, 0, 1, 0),
    return_tmb_report = FALSE,
    ...) {

  assert_that(inherits(object, "sdmTMB"))
  assert_that(length(nsim) == 1L)
  assert_that(is.numeric(nsim))
  assert_that(nsim >= 1)
  assert_that(length(nproj) == 1L)
  assert_that(is.numeric(nproj))
  assert_that(nproj >= 1)
  assert_that(is.data.frame(newdata))
  assert_that(length(sims_var) == 1L)
  assert_that(is.logical(return_tmb_report))
  assert_that(is.numeric(model))
  assert_that(length(model) == 1L)
  assert_that(as.integer(model) %in% c(1L, 2L))
  assert_that(is.logical(silent))
  assert_that(length(sim_re) == 6L)
  assert_that(as.integer(sim_re) %in% c(0L, 1L))

  ## extend time keeping elements of sdmTMB object
  max_year_i <- max(object$time_lu$year_i)
  new_year_i <- seq(max_year_i + 1, max_year_i + nproj)
  max_time <- max(object$extra_time)
  new_time_from_data <- seq(max_time + 1, max_time + nproj)
  object$extra_time <- c(object$extra_time, new_time_from_data)
  object$time_lu$sim_projected <- FALSE
  object$time_lu <- rbind(
    object$time_lu,
    data.frame(
      year_i = new_year_i, time_from_data = new_time_from_data,
      extra_time = TRUE, sim_projected = TRUE
    )
  )

  ## generate prediction TMB data list
  p <- predict(object, newdata = newdata, return_tmb_data = TRUE, ...)

  ## move data elements over
  p <- move_proj_to_tmbdat(p, object, newdata)

  ## extend time
  p$n_t <- nrow(object$time_lu)
  p$simulate_t <- as.integer(object$time_lu$sim_projected)
  ## sim random effects? order: omega, epsilon, zeta, IID, RW, smoothers
  p$sim_re <- sim_re

  ## parameters: add zeros as needed to all time-based parameters
  pars <- get_pars(object)
  n_s <- dim(pars$epsilon_st)[1]
  nt_new <- 1L
  pars$epsilon_st <- abind::abind(pars$epsilon_st, array(0, c(n_s, nproj, 1)), along = 2)
  pars$b_rw_t <- abind::abind(pars$b_rw_t, array(0, c(nproj, 1, 1)), along = 1)

  ## extend mapping as need
  ## FIXME: drop dims to 0 to avoid mapping complexity
  map <- object$tmb_map
  if ("b_rw_t" %in% names(map)) {
    map$b_rw_t <- c(map$b_rw_t, factor(rep(NA, nproj)))
  }
  if ("epsilon_st" %in% names(map)) {
    cli::cli_abort("This function hasn't been set up to work with maps in epsilon_st yet.")
  }

  ## rebuild TMB object
  obj <- TMB::MakeADFun(
    data = p,
    profile = object$control$profile,
    parameters = pars,
    map = map,
    random = object$tmb_random,
    DLL = "sdmTMB",
    silent = TRUE
  )

  ## do simulations
  if (!silent) cli::cli_progress_bar("Simulating projections", total = nsim)
  ret <- list()
  for (i in seq_len(nsim)) {
    if (!silent) cli::cli_progress_update()
    ret[[i]] <- obj$simulate()
  }
  cli::cli_process_done()
  if (return_tmb_report) {
    return(ret)
  }
  ret <- lapply(ret, \(x) x[[sims_var]][, model])
  do.call(cbind, ret)
}

move_proj_to_tmbdat <- function(x, object, newdata) {
  x$do_predict <- 0L
  x$year_i <- x$proj_year
  ## FIXME:
  ## x$A_st <- x$proj_mesh
  ## .cpp uses unique locations in projection but not in fitting:
  xy_cols <- object$spde$xy_cols
  proj_mesh <- fmesher::fm_basis(object$spde$mesh, loc = as.matrix(newdata[, xy_cols, drop = FALSE]))
  x$A_st <- proj_mesh
  x$A_spatial_index <- seq_len(dim(proj_mesh)[1]) - 1L
  x$X_threshold <- x$proj_X_threshold
  x$X_ij <- x$proj_X_ij
  x$X_rw_ik <- x$proj_X_rw_ik # FIXME: matches proj?
  x$z_i <- x$proj_z_i
  x$Zs <- x$proj_Zs
  x$Xs <- x$proj_Xs
  x$RE_indexes <- x$proj_RE_indexes
  x$offset_i <- x$proj_offset_i
  x$y_i <- matrix(NA, ncol = 1, nrow = nrow(x$proj_X_ij[[1]])) # fake
  x$weights_i <- rep(1, nrow(x$y_i)) # fake
  x$area_i <- rep(1, nrow(x$y_i)) # fak
  unique_size <- unique(x$size)
  if (length(unique_size) != 1L) {
    cli::cli_abort("This function hasn't been set up to work with binomial size specified yet.")
  }
  x$size <- rep(1, nrow(x$y_i)) # FIXME:
  x
}

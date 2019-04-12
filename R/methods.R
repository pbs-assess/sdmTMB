# dat <- sim(time_steps = 10, plot = TRUE, initial_betas = c(0.1, 0.2, -0.1), sigma_V = c(0,  0.2, 0.2))
# spde <- make_spde(x = dat$x, y = dat$y, n_knots = 50)
# plot_spde(spde)
# m <- sdmTMB(
#   data = dat, formula = observed ~ cov1, time = "time", include_spatial = T,
#   time_varying = ~ 0 + cov2 + cov2,
#   family = gaussian(link = "identity"), spde = spde
# )

#' @export
#' @import methods
print.sdmTMB <- function(x, ...) {
  # if (isTRUE(x$args$spatial_only)) {
  #   title <- "Spatial model fit by ML ['sdmTMB']\n"
  # } else {
  title <- "Spatiotemporal/spatial model fit by ML ['sdmTMB']\n"
  # }
  formula <- paste0("Formula: ", deparse(x$formula), "\n")
  # data <- paste0("Data: ", x$args$data, "\n")
  family <- paste0("Family: ", paste0(x$family$family, "(link = '", x$family$link, "')"), "\n")
  criterion <- paste0("ML criterion at convergence: ", round(x$model$objective, 3), "\n")
  fe_names <- colnames(model.matrix(x$formula, x$data))

  r <- x$tmb_obj$report()
  pars <- x$model$par
  b_j <- round(unname(pars[grep("b_j", names(pars))]), 2L)

  phi <- round(exp(as.list(pars)$ln_phi), 2L)
  range <- round(r$range, 2L)

  pre <- "Spatial SD (sigma_O): "
  if (!is.null(r$sigma_O)) {
    sigma_O <- paste0(pre, round(r$sigma_O, 2L), "\n")
  } else {
    sigma_O <- paste0(pre, "not estimated\n")
  }

  pre <- "Spatiotemporal SD (sigma_E): "
  if (!is.null(r$sigma_E)) {
    sigma_E <- paste0(pre, round(r$sigma_E, 2L), "\n")
  } else {
    sigma_E <- paste0(pre, "not estimated\n")
  }

  pre <- "Spatiotemporal AR1 correlation (rho): "
  if (!is.null(r$rho) && r$rho != 0L) {
    rho <- paste0(pre, round(r$rho, 2L), "\n")
  } else {
    rho <- ""
  }

  sr <- x$sd_report
  sr_se <- summary(sr)[,"Std. Error"]
  sr_est <- summary(sr)[,"Std. Error"]

  b_j_se <- unname(round(sr_se[grep("b_j", names(sr_se))], 2L))
  b_j <- unname(round(sr_se[grep("b_j", names(sr_est))], 2L))

  mm <- cbind(b_j, b_j_se)
  colnames(mm) <- c("coef.est", "coef.se")
  row.names(mm) <- fe_names
  mm

  cat(title,
    # formula,
    # data,
    family,
    sep = "")

  print(mm)

  cat("\n",
    paste0("Matern range parameter: ", range, "\n"),
    paste0("Dispersion parameter: ", phi, "\n"),
    sigma_O,
    sigma_E,
    rho,
    criterion,
    sep = ""
  )
}

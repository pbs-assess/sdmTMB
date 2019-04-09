# dat <- sim(time_steps = 5, plot = TRUE)
# spde <- make_spde(x = dat$x, y = dat$y, n_knots = 50)
# plot_spde(spde)
# m <- sdmTMB(
#   data = dat, formula = observed ~ 1, time = "time", include_spatial = F,
#   family = gaussian(link = "identity"), spde = spde
# )
#
# title <- "Spatiotemporal model fit by ML ['sdmTMB']\n"
#
# formula <- paste0("Formula: ", deparse(m$formula), "\n")
#
# data <- paste0("Data: ", m$args$data, "\n")
#
# family <- paste0("Family: ", deparse(m$args$family), "\n")
#
# criterion <- paste0("ML criterion at convergence: ", round(m$model$objective, 3), "\n")
#
# fe_names <- colnames(model.matrix(m$formula, m$data))
#
# r <- m$tmb_obj$report()
# pars <- as.list(m$model$par)
# b_j <- setNames(round(pars[[grep("b_j", names(pars))]], 3L), fe_names)
# phi <- setNames(round(exp(pars$ln_phi), 3L), "phi")
# range <- setNames(round(r$range, 3L), "range")
#
# if (!is.null(r$sigma_O)) {
#   sigma_O <- setNames(round(r$sigma_O, 3L), "range")
# } else {
#   sigma_O <- ""
# }
#
# if (!is.null(r$sigma_E)) {
#   sigma_E <- setNames(round(r$sigma_E, 3L), "sigma_E")
# } else {
#   sigma_E <- ""
# }
#
# if (!is.null(r$rho)) {
#   rho <- setNames(round(r$rho, 3L), "rho")
# } else {
#   rho <- ""
# }
#
# cat(title,
#   formula,
#   data,
#   family,
#   criterion
# )
#
# sr <- m$sd_report
# sr_se <- summary(sr)[,"Std. Error"]
#
# b_j_se <- setNames(round(sr_se[[grep("b_j", names(sr_se))]], 3L), fe_names)
#
# # Error terms:

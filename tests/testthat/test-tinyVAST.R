# if (FALSE) {
# if (requireNamespace("tinyVAST", quietly = TRUE)) {
#   library("tinyVAST", warn.conflicts = FALSE)
#   TOL <- 1e-5
#
#   get_sdmTMB_pars <- function(x) {
#     c(
#       as.list(x$sd_report, "Estimate"),
#       as.list(x$sd_report, "Estimate", report = TRUE)
#     )
#   }
#
#   test_that("tinyVAST/sdmTMB Tweedie spatiotemporal IID models and index area integration match", {
#     skip_on_cran()
#     mesh <- make_mesh(pcod, c("X", "Y"), cutoff = 18)
#     system.time({
#       fit_sd <- sdmTMB(
#         density ~ 0 + as.factor(year),
#         data = pcod,
#         mesh = mesh,
#         family = tweedie(),
#         time = "year",
#         control = sdmTMBcontrol(multiphase = FALSE)
#       )
#     })
#     ps <- get_sdmTMB_pars(fit_sd)
#
#     system.time({
#       fit_tv <- tinyVAST(
#         density ~ 0 + factor(year),
#         dsem = "",
#         sem = "",
#         data = pcod,
#         family = tweedie(),
#         time_column = "year",
#         space_columns = c("X", "Y"),
#         spatial_graph = mesh$mesh,
#         control = tinyVASTcontrol(newton_loops = 1)
#       )
#     })
#     pt <- as.list(fit_tv$sdrep, "Estimate")
#
#     expect_equal(pt$alpha_j, ps$b_j, tolerance = TOL)
#     expect_equal(pt$log_sigma[2], ps$thetaf, tolerance = TOL)
#     expect_equal(pt$log_sigma[1], ps$ln_phi, tolerance = TOL)
#     expect_equal(pt$beta_z, ps$sigma_E[1, 1], tolerance = TOL)
#     expect_equal(pt$theta_z, ps$sigma_O[1, 1], tolerance = TOL)
#     expect_equal(pt$log_kappa, ps$ln_kappa[1, 1], tolerance = TOL)
#
#     g <- replicate_df(qcs_grid, "year", unique(pcod$year))
#     p <- predict(fit_sd, newdata = g, return_tmb_object = TRUE)
#
#     system.time({
#       is <- get_index(p, bias_correct = TRUE)
#     })
#
#     # system.time({
#     # is2 <- lapply(unique(g$year), \(x) {
#     #   pp <- predict(fit_sd, newdata = subset(g, year == x), return_tmb_object = TRUE)
#     #   get_index(pp, bias_correct = TRUE)
#     # })
#     # is2 <- do.call(rbind, is2)
#     # })
#
#     system.time({
#       g$var <- "response"
#       it <- lapply(unique(g$year), \(x)
#       integrate_output(fit_tv, newdata = subset(g, year == x), apply.epsilon = TRUE))
#       it <- do.call(rbind, it) |> as.data.frame()
#     })
#
#     expect_equal(it$`Est. (bias.correct)`, is$est, tolerance = TOL)
#   })
#
#   fit_sd <- sdmTMB(
#     density ~ 0 + as.factor(year),
#     data = pcod,
#     family = delta_lognormal(type = "poisson-link"),
#     spatial = "off",
#     spatiotemporal = "off"
#   )
#   ps <- get_pars(fit_sd)
#
#   suppressWarnings({
#     fit_tv <- tinyVAST(
#       density ~ 0 + factor(year),
#       delta_options = list(delta_formula = ~ 0 + factor(year)),
#       data = pcod,
#       family = delta_lognormal(type = "poisson-link"),
#       control = tinyVASTcontrol(newton_loops = 1)
#     )
#   })
#   pt <- as.list(fit_tv$sdrep, "Estimate")
#
#   compare_models_delta <- function(mtv, msd) {
#     ps <- get_pars(msd)
#     pt <- as.list(mtv$sdrep, "Estimate")
#     expect_equal(logLik(msd), logLik(mtv), ignore_attr = TRUE)
#     expect_equal(pt$alpha_j, ps$b_j, tolerance = TOL)
#     expect_equal(pt$alpha2_j, ps$b_j2, tolerance = TOL)
#   }
#
#   test_that("tinyVAST/sdmTMB delta_lognormal(type = 'poisson-link') models match", {
#     skip_on_cran()
#     compare_models_delta(fit_tv, fit_sd)
#     expect_equal(pt$log_sigma, ps$ln_phi[2], tolerance = TOL)
#   })
#
#   test_that("tinyVAST/sdmTMB delta_gamma() models match", {
#     skip_on_cran()
#     fit_sd <- update(fit_sd, family = delta_gamma())
#     fit_tv <- update(fit_tv, family = delta_gamma())
#     compare_models_delta(fit_tv, fit_sd)
#     ps <- get_pars(fit_sd)
#     pt <- as.list(fit_tv$sdrep, "Estimate")
#     expect_equal(1 / (exp(pt$log_sigma))^2, exp(ps$ln_phi[2]), tolerance = TOL)
#   })
#
#   test_that("tinyVAST/sdmTMB delta_gamma(type = 'poisson-link')) models match", {
#     skip_on_cran()
#     fit_sd <- update(fit_sd, family = delta_gamma(type = "poisson-link"))
#     fit_tv <- update(fit_tv, family = delta_gamma(type = "poisson-link"))
#     compare_models_delta(fit_tv, fit_sd)
#     ps <- get_pars(fit_sd)
#     pt <- as.list(fit_tv$sdrep, "Estimate")
#     expect_equal(1 / (exp(pt$log_sigma))^2, exp(ps$ln_phi[2]), tolerance = TOL)
#   })
# }
# }
#

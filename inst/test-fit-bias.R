# context("Is there any bias in basic fits?")
#
# SEED <- 1
# set.seed(SEED)
# x <- stats::runif(100, -1, 1)
# y <- stats::runif(100, -1, 1)
#
# # test_that("sdmTMB model fit with a covariate beta", {
# initial_betas <- 0.5
# kappa <- 6 # decay of spatial correlation (smaller = slower decay)
# sigma_O <- 0.3 # SD of spatial process
# sigma_E <- 0.3 # SD of spatial process
# phi <- 0.05 # observation error
#
# library(future)
# plan(multisession, workers = 8L)
#
# library(sdmTMB)
# # out <- furrr::future_map(seq(1, 8*2), function(i) {
# out <- purrr::map(seq(1, 8*2), function(i) {
#   # s <- sim(
#   #   x = x, y = y,
#   #   initial_betas = initial_betas, time_steps = 10L,
#   #   phi = phi, kappa = kappa, sigma_O = sigma_O, sigma_E = sigma_E,
#   #   seed = SEED * i
#   # )
#   set.seed(SEED * i)
#   x <- stats::runif(100, -1, 1)
#   y <- stats::runif(100, -1, 1)
#   s <- sim(
#     x = x, y = y,
#     initial_betas = initial_betas, time_steps = 1L,
#     phi = phi, kappa = kappa, sigma_O = sigma_O,
#     seed = SEED * i
#   )
#
#   # coords <- cbind(x, y)
#   # new <- list(n = Inf)
#   # cutoff <- 0.1
#   # counter <- 0
#   # percent <- 0.06
#   # increment <- .005
#   # numknots <- 100
#   # while (new$n > numknots * (1 + percent) || new$n < numknots * (1 - percent)) {
#   #   cat(counter, "|| cutoff =", sprintf("%.2f", cutoff), "\n")
#   #   cutoff <- ifelse(new$n > numknots * (1 + percent),
#   #     cutoff + increment, cutoff - increment)
#   #   counter <- counter + 1
#   #   new <- INLA::inla.mesh.create(coords, refine = TRUE, cutoff = cutoff)
#   #   if (counter > 100 || cutoff < 0.001) break
#   # }
#   # print(cutoff)
#   # print(new$n)
#   #
#   # new <- INLA::inla.mesh.create(coords, refine = TRUE, cutoff = 0.001)
#   # spde <- make_spde(s$x, s$y, mesh = new)
#   # plot_spde(spde);points(coords, col = "blue")
#   #
#   # spde <- make_spde(x = s$x, y = s$y, n_knots = 70)
#   # plot_spde(spde);points(coords, col = "blue")
#   #
#   loc_xy <- cbind(s$x, s$y)
#   bnd <- INLA::inla.nonconvex.hull(as.matrix(loc_xy), convex = -0.05)
#   mesh <- INLA::inla.mesh.2d(
#     boundary = bnd,
#     max.edge = c(20, 20),
#     offset = -0.05,
#     cutoff = c(0.01, 0.1),
#     min.angle = 10
#   )
#   sp2 <- make_spde(s$x, s$y, mesh = mesh)
#   plot_spde(sp2);points(coords, col = "blue")
#   spde <- sp2
#   # m <- tryCatch({sdmTMB(data = s, formula = observed ~ 0 + cov1,
#     # time = "time", spde = spde)}, error = function(e) NA)
#   m <- tryCatch({sdmTMB(data = s, formula = observed ~ 0 + cov1,
#     time = "time", spatial_only = TRUE, spde = spde)}, error = function(e) NA)
#   if (identical(m, NA)) return(NA)
#   p <- as.list(m$model$par)
#   r <- m$tmb_obj$report()
#   data.frame(sigma_O = r$sigma_O,
#     # sigma_E = r$sigma_E,
#     b = p$b_j,
#     kappa = exp(p$ln_kappa), phi = exp(p$ln_phi))
# })
#
# library(dplyr)
# d <- dplyr::bind_rows(out)
# true <- tibble(variable = c("sigma_O", "sigma_E", "b", "kappa", "phi"),
#   true_value = c(sigma_O, sigma_E, initial_betas, kappa, phi))
# library(ggplot2)
# reshape2::melt(d) %>% left_join(true) %>% ggplot(aes(value)) +
#   facet_wrap(vars(variable), scales = "free") +
#   geom_histogram(bins = 20) +
#   geom_vline(aes(xintercept = true_value), colour = "red")
#
#
# plan(sequential)
#
# # expect_equal(m$model$convergence, 0L)
# # expect_equal((p$b_j - initial_betas)^2, 0, tol = 0.05)
# # expect_equal((exp(p$ln_phi) - phi)^2, 0, tol = 0.05)
# # expect_equal((r$sigma_O - sigma_O)^2, 0, tol = 0.05)
# # expect_equal((r$sigma_E - sigma_E)^2, 0, tol = 0.05)
# # expect_equal(exp(p$ln_kappa), kappa, tol = 1.1)
# # p <- predict(m)
# # r <- residuals(m)
# # expect_equal(mean((p$est - s$observed)^2), 0, tol = 0.02)
# # })

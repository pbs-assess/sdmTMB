# basic model fitting and prediction tests

test_that("sdmTMB model estimates standard curves", {
  local_edition(2)

  set.seed(123)
  known_conc_ul <- c(1e05, 1e04, 1e03, 1e02, 1e01, 1, 5, 2.322919, 2.050185, 1.289870e-01, 4.922461e-02, 1.633000e+02, 4.598000, 4.799000, 2.897600, 4.387000, 4.185000, 6.369000)

  plates = paste0(sample(letters, 50, replace=T), sample(1:100, size=50))
  standards <- expand.grid(known_conc_ul = known_conc_ul, replicate = 1:3, plate = plates)
  standards$plate <- as.factor(standards$plate)

  Sigma <- matrix(c(0.26, 0.081, 0.081, 0.026), 2, 2)
  logistic_coefs <- mvtnorm::rmvnorm(n = length(plates), c(2,0.1), Sigma)
  Sigma <- matrix(c(0.65, -0.01, -0.01, 0.00145), 2, 2)
  gauss_coefs <- mvtnorm::rmvnorm(n = length(plates), c(30,-1.5), Sigma)
  true_coefs <- as.data.frame(cbind(logistic_coefs, gauss_coefs))
  names(true_coefs) = c("std_xi_2_true", "std_xi_3_true", "std_xi_0_true", "std_xi_1_true")

  phi <- 0.01
  standards$Ct <- gauss_coefs[match(standards$plate, levels(standards$plate)),1] + gauss_coefs[match(standards$plate, levels(standards$plate)),2] * log(standards$known_conc_ul) + rnorm(nrow(standards), 0, phi)

  # presence - absence model
  standards$p <- plogis(logistic_coefs[match(standards$plate, levels(standards$plate)),1] + logistic_coefs[match(standards$plate, levels(standards$plate)),2] * log(standards$known_conc_ul))
  standards$Ct <- standards$Ct * ifelse(runif(nrow(standards)) < standards$p, 1, 0)

  predictor_dat <- data.frame(
    X = runif(500, 0, 500), Y = runif(500, 0, 500), year = 1
  )
  mesh <- make_mesh(predictor_dat, xy_cols = c("X", "Y"), cutoff = 20)

  sim_dat <- sdmTMB_simulate(
    formula = ~ 1,
    data = predictor_dat,
    time = "year",
    mesh = mesh,
    family = gaussian(),
    range = 0.5,
    sigma_E = 0,
    phi = 0.01,
    sigma_O = 0.2,
    seed = 123,
    B = 3 # B0 = intercept, B1 = a1 slope
  )

  sim_dat$plate <- sample(unique(standards$plate), nrow(sim_dat), replace=T)

  sim_dat$Ct <- gauss_coefs[match(sim_dat$plate, levels(standards$plate)),1] + gauss_coefs[match(sim_dat$plate, levels(standards$plate)),2] * sim_dat$eta + rnorm(nrow(sim_dat), 0, phi)

  # presence - absence model
  sim_dat$p <- plogis(logistic_coefs[match(sim_dat$plate, levels(standards$plate)),1] + logistic_coefs[match(sim_dat$plate, levels(standards$plate)),2] * sim_dat$eta)
  sim_dat$Ct <- sim_dat$Ct * ifelse(runif(nrow(sim_dat)) < sim_dat$p, 1, 0)

  mesh <- make_mesh(sim_dat, xy_cols = c("X", "Y"), cutoff = 25)
  fit <- sdmTMB(Ct ~ 1,
                mesh = mesh,
                spatial = "on",
                spatiotemporal = "off",
                control = sdmTMBcontrol(stdcurve_df = standards),
                data=sim_dat)

  expect_equal(as.numeric(tidy(fit)[,"estimate"]), 2.9997, tolerance=0.0001)

  r <- fit$sd_report$par.random

  df <- data.frame(id = fit$plates$id,
                   plate = fit$plates$plate,
                   std_xi_0 = r[grep("std_xi_0", names(r))],
                   std_xi_1 = r[grep("std_xi_1", names(r))],
                   std_xi_2 = r[grep("std_xi_2", names(r))],
                   std_xi_3 = r[grep("std_xi_3", names(r))]
  )

  # bring in true coefficients
  df <- cbind(df, true_coefs)

  expect_equal(cor(df$std_xi_0, df$std_xi_0_true), 1, tolerance=0.0001)
  expect_equal(cor(df$std_xi_1, df$std_xi_1_true), 1, tolerance=0.0001)
  expect_gt(cor(df$std_xi_2, df$std_xi_2_true), 0.6)
  expect_gt(cor(df$std_xi_3, df$std_xi_3_true), 0.6)

  # remove plates from df
  fit <- try(sdmTMB(Ct ~ 1,
                mesh = mesh,
                spatial = "on",
                spatiotemporal = "off",
                control = sdmTMBcontrol(stdcurve_df = standards),
                data=dplyr::select(sim_dat,-plate)), silent=TRUE)
  expect_equal(class(fit), "try-error")

  fit <- try(sdmTMB(Ct ~ 1,
                    mesh = mesh,
                    spatial = "on",
                    spatiotemporal = "off",
                    control = sdmTMBcontrol(stdcurve_df = dplyr::select(standards,-plate)),
                    data=dplyr::select(sim_dat)), silent=TRUE)
  expect_equal(class(fit), "try-error")

  fit <- try(sdmTMB(Ct ~ 1,
                    mesh = mesh,
                    spatial = "on",
                    spatiotemporal = "off",
                    control = sdmTMBcontrol(stdcurve_df = dplyr::select(standards,-Ct)),
                    data=dplyr::select(sim_dat)), silent=TRUE)
  expect_equal(class(fit), "try-error")

  fit <- try(sdmTMB(Ct ~ 1,
                    mesh = mesh,
                    spatial = "on",
                    spatiotemporal = "off",
                    control = sdmTMBcontrol(stdcurve_df = dplyr::select(standards,-known_conc_ul)),
                    data=dplyr::select(sim_dat)), silent=TRUE)
  expect_equal(class(fit), "try-error")
})


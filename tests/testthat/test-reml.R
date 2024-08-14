test_that("REML works", {
  skip_on_cran()

  set.seed(1)

  pcod$year_f <- as.factor(pcod$year)

  fit_glmm <- glmmTMB::glmmTMB(present ~ 1 + (1 | year_f),
    data = pcod,
    family = binomial()
  )

  fit_glmm_reml <- glmmTMB::glmmTMB(present ~ 1 + (1 | year_f),
    data = pcod,
    REML = TRUE,
    family = binomial()
  )

  fit_sdm <- sdmTMB(present ~ 1 + (1 | year_f),
    data = pcod, spatial = "off",
    family = binomial()
  )

  # warning thrown if glmmTMB argument mistakenly used
  expect_error({
    fit_sdm_wrong <- sdmTMB(present ~ 1 + (1 | year_f),
      data = pcod, spatial = "off",
      REML = TRUE,
      family = binomial()
    )
  })

  fit_sdm_reml <- sdmTMB(present ~ 1 + (1 | year_f),
    data = pcod, spatial = "off",
    reml = TRUE,
    family = binomial()
  )


  AIC(fit_glmm, fit_sdm, fit_glmm_reml, fit_sdm_reml)

  # check that the AIC calculations match
  expect_equal(AIC(fit_glmm), AIC(fit_sdm), tolerance = 1e-8)
  expect_equal(AIC(fit_glmm_reml), AIC(fit_sdm_reml), tolerance = 1e-8)

  # check that the global intercepts match for ML
  b1 <- summary(fit_glmm)$coefficients$cond[1]
  b2 <- tidy(fit_sdm)$estimate[1]

  expect_equal(b1, b2, tolerance = 1e-5)

  # check that the global intercepts still match for REML
  b1_reml <- summary(fit_glmm_reml)$coefficients$cond[1]
  b2_reml <- tidy(fit_sdm_reml)$estimate[1]

  expect_equal(b1_reml, b2_reml, tolerance = 1e-5)

  # check that the random effect variances match for ML
  b1_re <- sqrt(summary(fit_glmm)$varcor$cond[[1]])[1]
  b2_re <- tidy(fit_sdm, effects = "ran_pars")$estimate[1]

  expect_equal(b1_re, b2_re, tolerance = 1e-5)

  # check that the random effect variances still match for REML

  b1_reml_re <- sqrt(summary(fit_glmm_reml)$varcor$cond[[1]])[1]
  b2_reml_re <- tidy(fit_sdm_reml, effects = "ran_pars")$estimate[1]

  expect_equal(b1_reml_re, b2_reml_re, tolerance = 1e-5)
})

test_that("REML works for delta models", {
  skip_on_cran()

  set.seed(1)

  pcod$year_f <- as.factor(pcod$year)
  pcod_pos <- subset(pcod, density > 0)


  fit_glmm <- glmmTMB::glmmTMB(density ~ 1 + (1 | year_f),
    data = pcod_pos,
    family = Gamma("log")
  )

  fit_sdm <- sdmTMB(density ~ 1 + (1 | year_f),
    data = pcod_pos,
    spatial = "off",
    family = Gamma("log")
  )

  fit_glmm_reml <- glmmTMB::glmmTMB(density ~ 1 + (1 | year_f),
    data = pcod_pos,
    REML = TRUE,
    family = Gamma("log")
  )

  fit_sdm_reml <- sdmTMB(density ~ 1 + (1 | year_f),
    data = pcod_pos,
    spatial = "off",
    reml = TRUE,
    family = Gamma("log")
  )

  AIC(fit_glmm, fit_sdm, fit_glmm_reml, fit_sdm_reml)

  # check that the AIC calculations match
  expect_equal(AIC(fit_glmm), AIC(fit_sdm), tolerance = 1e-8)
  expect_equal(AIC(fit_glmm_reml), AIC(fit_sdm_reml), tolerance = 1e-8)

  # check that the global intercepts match for ML
  b1 <- summary(fit_glmm)$coefficients$cond[1]
  b2 <- tidy(fit_sdm)$estimate[1]

  expect_equal(b2, b1, tolerance = 1e-6)

  # check that the global intercepts still match for REML
  b1_reml <- summary(fit_glmm_reml)$coefficients$cond[1]
  b2_reml <- tidy(fit_sdm_reml)$estimate[1]

  expect_equal(b2_reml, b1_reml, tolerance = 1e-6)

  # REML is doing something more than just difference between model runs
  expect_false(((b2_reml - b2) > 1e-4))

  # now try delta models
  fit_dg_reml <- sdmTMB(density ~ 1 + (1 | year_f),
    data = pcod, spatial = "off",
    reml = TRUE,
    family = delta_gamma()
  )

  # global intercepts match delta model 2 with REML
  b3_reml <- tidy(fit_dg_reml, model = 2)$estimate[1]
  expect_equal(b3_reml, b2_reml, tolerance = 1e-6)

  # sigma_G match delta model 2 with REML
  b2_reml_re <- tidy(fit_sdm_reml, effects = "ran_pars")
  b3_reml_re <- tidy(fit_dg_reml, effects = "ran_pars", model = 2)

  expect_equal(b3_reml_re$estimate[b3_reml_re$term=="sigma_G"],
               b2_reml_re$estimate[b2_reml_re$term=="sigma_G"], tolerance = 1e-6)

  # REML is still doing something with delta models with spatial fields
  mesh <- make_mesh(pcod, c("X", "Y"), cutoff = 10)
  fit_dg2 <- sdmTMB(density ~ 1 + (1 | year_f),
    data = pcod,
    spatial = "on", mesh = mesh,
    family = delta_gamma()
  )

  fit_dg2_reml <- sdmTMB(density ~ 1 + (1 | year_f),
    data = pcod,
    spatial = "on", mesh = mesh,
    reml = TRUE,
    family = delta_gamma()
  )

  b4 <- tidy(fit_dg2, model = 2)$estimate[1]
  b4_reml <- tidy(fit_dg2_reml, model = 2)$estimate[1]

  # effect of reml on fixed effects in delta model
  expect_false((abs(b4_reml - b4) < 0.01))

  b4_re <- tidy(fit_dg2, effects = "ran_pars", conf.int = T, model = 2)#$estimate[4]
  b4_reml_re <- tidy(fit_dg2_reml, effects = "ran_pars", conf.int = T, model = 2)#$estimate[4]

  # effect of reml on sigma_G in delta model
  expect_false((abs(b4_reml_re$estimate[b4_reml_re$term=="sigma_G"] - b4_re$estimate[b4_re$term=="sigma_G"]) < 0.01))
})

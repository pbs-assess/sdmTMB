test_that("get_index(), get_index_sims(), and get_cog() work", {
  skip_on_ci()
  skip_on_cran()
  skip_if_not_installed("INLA")
  pcod_spde <- make_mesh(pcod, c("X", "Y"), cutoff = 20)
  m <- sdmTMB(
    data = pcod,
    formula = density ~ 0 + as.factor(year),
    time = "year", spde = pcod_spde, family = tweedie(link = "log")
  )
  predictions <- predict(m, newdata = qcs_grid, return_tmb_object = TRUE)

  p <- predict(m, newdata = qcs_grid, return_tmb_object = FALSE)
  expect_error(get_index(p), regexp = "return_tmb_object")

  ind <- get_index(predictions, bias_correct = FALSE)

  expect_equal(class(ind), "data.frame")
  expect_equal(ind$est,
    c(209277.00119385, 316826.468229804, 355481.950149426, 91540.9936646985,
    150187.860588992, 274230.962985684, 274541.848909578, 326960.546698788,
    151502.369101489))
  expect_equal(ind$lwr,
    c(151512.996010835, 244648.356445529, 274276.862571144, 66923.2486410552,
    109088.999957795, 212808.468150465, 204421.906823893, 241466.295497819,
    110227.656055564))
  expect_equal(ind$upr,
    c(289063.40962038, 410299.142938737, 460729.409318153, 125214.386499159,
    206770.558690844, 353381.713207407, 368714.038401094, 442725.138418069,
    208232.38617896))
  set.seed(1)
  pred_sim <- predict(m, sims = 600L)
  ind_sim <- get_index_sims(pred_sim)
  expect_gt(cor(ind_sim$est, ind$est), 0.95)
  expect_gt(cor(ind_sim$lwr, ind$lwr), 0.95)
  expect_gt(cor(ind_sim$upr, ind$upr), 0.95)
  # sims mimics bias corrected one, which would be higher:
  expect_gt(mean(ind$est), mean(ind_sim$est))

  cog <- get_cog(predictions, bias_correct = FALSE)
  expect_equal(class(cog), "data.frame")
  expect_equal(cog$est,
    c(462.086501154086, 481.335275610853, 471.671456899259, 481.90130273488,
    485.463933915529, 469.799589078202, 475.98239867458, 457.43341815583,
    463.397343260075, 5756.9908955009, 5729.23474156994, 5760.98573470148,
    5735.0944675103, 5726.17650067701, 5745.24969598697, 5744.28285295928,
    5757.13546800737, 5755.77780115423))
})


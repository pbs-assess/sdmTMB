# comparing VAST and sdmTMB
#To do:
  #Distributions: Tweedie, delta-gamma
    #Poisson-link delta
  #Random walk
  #AR1
  #IID
  #Spatially varying coefficient

test_that("VAST test", {
  #VAST setup
  install.packages("devtools")
  library(devtools)
  library(VAST)
  #sdmTMB setup
  install.packages("sdmTMB")
  library(sdmTMB)
  #sdmTMB model
  spde <- make_mesh(loc, c("x", "y"), cutoff = 10)
  fit <- sdmTMB()

  #VAST model
  expect_equal(fit$AIC, fit2$tolerance=)
  expect_equal()  #annual index
  })


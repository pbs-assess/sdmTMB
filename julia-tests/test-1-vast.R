# comparing VAST and sdmTMB
#To do:
  #Distributions: Tweedie, delta-gamma, delta-lognormal
  #Spatiotemporal structures: IID, AR1, RW
  #Parameter estimates, predictions, COG, index values
  #Bias correction

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
  fit <- sdmTMB(
    density ~ depth,
    mesh= mesh,
    family=tweedie(link="log"),
    spatial="on",
    spatiotemporal="on"
  )
  #VAST model
  fit2  <-
  expect_equal(fit$AIC, fit2$AIC, tolerance=)
  expect_equal()  #annual index
  })


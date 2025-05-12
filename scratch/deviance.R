# library(sdmTMB)
# mesh <- make_mesh(pcod_2011, c("X", "Y"), cutoff = 20)
# fit <- sdmTMB(
#   density ~ depth_scaled,
#   data = pcod_2011, mesh = mesh,
#   family = tweedie(link = "log")
# )
# fit
# r <- fit$tmb_obj$report()

x <- rnorm(10)
y <- rpois(10, exp(x))
d <- data.frame(x = x, y = y)
m <- sdmTMB(y ~ 1 + x, family = poisson(), spatial = "off", data = d)
mglm <- glm(y ~ 1 + x, family = poisson(), data = d)
r1 <- unname(residuals(mglm))
rep <- m$tmb_obj$report()
r2 <- rep$devresid[, 1]
expect_equal(r1, r2, tolerance = 0.01)
deviance(mglm)
sum(r2^2)

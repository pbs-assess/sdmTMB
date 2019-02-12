d <- pcod
pcod_spde <- make_spde(d$X, d$Y, n_knots = 80) # only 100 knots for example speed
plot_spde(pcod_spde)

# Tweedie:
m_ar1 <- sdmTMB(
  d, density ~ 0 + as.factor(year) + depth_scaled + depth_scaled2,
  time = "year", spde = pcod_spde, family = tweedie(link = "log"),
  silent = FALSE, ar1_fields = TRUE)

m <- sdmTMB(
  d, density ~ 0 + as.factor(year) + depth_scaled + depth_scaled2,
  time = "year", spde = pcod_spde, family = tweedie(link = "log"),
  silent = FALSE, include_spatial = FALSE)

m_ar1$model$par

minus_one_to_one <- function(x) 2 * plogis(x) - 1
minus_one_to_one(m_ar1$model$par[["ar1_phi"]])


m$model$par

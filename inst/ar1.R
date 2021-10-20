d <- pcod
pcod_spde <- make_spde(d$X, d$Y, n_knots = 50)
plot_spde(pcod_spde)

# Tweedie:
m_ar1 <- sdmTMB(
  d, density ~ 0 + as.factor(year) + depth_scaled + depth_scaled2,
  time = "year", spde = pcod_spde, family = tweedie(link = "log"),
  silent = FALSE, spatiotemporal = "AR1", include_spatial = FALSE)

m <- sdmTMB(
  d, density ~ 0 + as.factor(year) + depth_scaled + depth_scaled2,
  time = "year", spde = pcod_spde, family = tweedie(link = "log"),
  silent = FALSE)

m$model$par
m_ar1$model$par

m_ar1$model$objective
m$model$objective

minus_one_to_one <- function(x) 2 * plogis(x) - 1
minus_one_to_one(m_ar1$model$par[["ar1_phi"]])

# ---------------------------------

set.seed(28331)
gp_sigma <- 0.2
sigma <- 0.1
df <- 1000
gp_theta <- 1.8
n_draws <- 20
nknots <- 40
phi <- 0.3
s <- glmmfields::sim_glmmfields(
  df = df, n_draws = n_draws, gp_theta = gp_theta,
  gp_sigma = gp_sigma, sd_obs = sigma, n_knots = nknots, phi = phi,
  n_data_points = 200
)
print(s$plot)
spde <- make_spde(s$dat$lon, s$dat$lat, n_knots = 60)
plot_spde(spde)

m_ar1_sim <- sdmTMB(
  s$dat, y ~ 0, time = "time", spde = spde,
  silent = FALSE, spatiotemporal = "AR1", spatial = FALSE)

m_ar1_sim$model$par
minus_one_to_one(m_ar1_sim$model$par[["ar1_phi"]])

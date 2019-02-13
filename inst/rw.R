d <- pcod
pcod_spde <- make_spde(d$X, d$Y, n_knots = 50)
plot_spde(pcod_spde)

# Tweedie:
m_rw1 <- sdmTMB(
  d, density ~ 0 + as.factor(year) + depth_scaled,
  time = "year", time_varying = "depth_scaled2",
  spde = pcod_spde, family = tweedie(link = "log"),
  silent = FALSE, ar1_fields = FALSE, include_spatial = TRUE)

# set.seed(28331)
#
# gp_sigma <- 0.2
# sigma <- 0.1
# df <- 1000
# gp_theta <- 1.8
# n_draws <- 12
# nknots <- 5
# year_sigma <- 0.5
# B <- vector(mode = "double", length = n_draws)
# B[1] <- 0
# for (i in 2:length(B)) {
#   B[i] <- B[i - 1] + rnorm(1, 0, year_sigma) # random walk
# }
#
# cov_vec = rnorm(n_draws*100,0,1)
# model_matrix = model.matrix(~a - 1 + cov + cov2,
#   data.frame(a = gl(n_draws, 100), cov=cov_vec, cov2=cov_vec^2))
#
# s <- glmmfields::sim_glmmfields(
#   df = df, n_draws = n_draws, gp_theta = gp_theta,
#   gp_sigma = gp_sigma, sd_obs = sigma, n_knots = nknots,
#   B = c(B, 3, -0.1),
#   X = model_matrix)
# s$dat$cov = model_matrix[,"cov"]
# s$dat$cov2 = model_matrix[,"cov2"]
#
# print(s$plot)
# spde <- make_spde(s$dat$lon, s$dat$lat, n_knots = 60)
# plot_spde(spde)
#
# m_rw <- sdmTMB(
#   s$dat, y ~ -1 + cov1 + cov2, time = "time", spde = spde,
#   silent = FALSE, ar1_fields = TRUE, include_spatial = FALSE)
#
# m_ar1_sim$model$par
# minus_one_to_one(m_ar1_sim$model$par[["ar1_phi"]])

devtools::load_all()
# library(ggplot2)

mesh <- make_mesh(dogfish, c("X", "Y"), cutoff = 25)
historical_years <- 2004:2022
to_project <- 10
future_years <- seq(max(historical_years) + 1, max(historical_years) + to_project)
all_years <- c(historical_years, future_years)
proj_grid <- replicate_df(wcvi_grid, "year", all_years)

fit2 <- sdmTMB(
  catch_weight ~ 1,
  time = "year",
  offset = log(dogfish$area_swept),
  # time_varying = ~ 1,
  # time_varying_type = "ar1",
  extra_time = historical_years, #< does *not* include projection years
  spatial = "on",
  spatiotemporal = list("off", "ar1"),
  data = dogfish,
  mesh = mesh,
  family = delta_gamma()
)

set.seed(1)
out <- project(fit2, newdata = proj_grid, nproj = to_project, nsim = 2, uncertainty = "joint")

est_se <- apply(out, 1, sd)
proj_grid$est_se_both <- est_se
ggplot(subset(proj_grid, year > 2021), aes(X, Y, fill = est_se)) +
  geom_raster() +
  facet_wrap(~year) +
  scale_fill_viridis_c() +
  coord_fixed() +
  ggtitle("Estimate SE (both)")

out <- project(fit2, newdata = proj_grid, nproj = to_project, nsim = 150, uncertainty = "none")
est_se <- apply(out, 1, sd)
proj_grid$est_se_none <- est_se
ggplot(subset(proj_grid, year > 2021), aes(X, Y, fill = est_se_none)) +
  geom_raster() +
  facet_wrap(~year) +
  scale_fill_viridis_c() +
  coord_fixed() +
  ggtitle("Estimate SE (both)")

out <- project(fit2, newdata = proj_grid, nproj = to_project, nsim = 150, uncertainty = "random")
est_se <- apply(out, 1, sd)
proj_grid$est_se_random <- est_se
ggplot(subset(proj_grid, year > 2021), aes(X, Y, fill = est_se_random)) +
  geom_raster() +
  facet_wrap(~year) +
  scale_fill_viridis_c() +
  coord_fixed() +
  ggtitle("Estimate SE (random)")

future_grid <- subset(proj_grid, year > 2022)
ggplot(future_grid, aes(est_se_random, est_se_none)) + geom_point() +
  geom_abline(intercept = 0, slope = 1, colour = "red", lwd = 1)
ggplot(future_grid, aes(est_se_both, est_se_none)) + geom_point() +
  geom_abline(intercept = 0, slope = 1, colour = "red", lwd = 1)
ggplot(future_grid, aes(est_se_both, est_se_random)) + geom_point() +
  geom_abline(intercept = 0, slope = 1, colour = "red", lwd = 1)


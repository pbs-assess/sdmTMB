## ----setup, include = FALSE, cache=FALSE--------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.width = 7,
  fig.asp = 0.618
)

## ----packages, message=FALSE, warning=TRUE------------------------------------
library(ggplot2)
library(dplyr)
library(sdmTMB)

## -----------------------------------------------------------------------------
pcod_spde <- make_spde(pcod$X, pcod$Y, n_knots = 100)
plot_spde(pcod_spde)

## -----------------------------------------------------------------------------
m1 <- sdmTMB(density ~ 1, data = pcod,
  spde = pcod_spde, family = tweedie(link = "log"),
  spatial_trend = TRUE, time = "year",
  spatial_only = TRUE)

## -----------------------------------------------------------------------------
m2 <- sdmTMB(density ~ 1, data = pcod,
  spde = pcod_spde, family = tweedie(link = "log"),
  spatial_trend = TRUE, time = "year",
  spatial_only = FALSE)

## -----------------------------------------------------------------------------
m3 <- sdmTMB(density ~ 1, data = pcod,
  spde = pcod_spde, family = tweedie(link = "log"),
  spatial_trend = TRUE, time = "year",
  spatial_only = FALSE, ar1_fields = TRUE)

## -----------------------------------------------------------------------------
d <- pcod
d$residuals1 <- residuals(m1)
d$residuals2 <- residuals(m2)
d$residuals3 <- residuals(m3)

qqnorm(d$residuals1);abline(a = 0, b = 1)
qqnorm(d$residuals2);abline(a = 0, b = 1)
qqnorm(d$residuals3);abline(a = 0, b = 1)

## -----------------------------------------------------------------------------
plot_map_point <- function(dat, column = "est") {
  ggplot(dat, aes_string("X", "Y", colour = column)) +
    geom_point() +
    facet_wrap(~year) +
    coord_fixed()
}

## -----------------------------------------------------------------------------
plot_map_point(d, "residuals1") + scale_color_gradient2()
plot_map_point(d, "residuals2") + scale_color_gradient2()
plot_map_point(d, "residuals3") + scale_color_gradient2()

## -----------------------------------------------------------------------------
sd1 <- as.data.frame(summary(TMB::sdreport(m1$tmb_obj)))
sd2 <- as.data.frame(summary(TMB::sdreport(m2$tmb_obj)))
sd3 <- as.data.frame(summary(TMB::sdreport(m3$tmb_obj)))

## -----------------------------------------------------------------------------
r1 <- m1$tmb_obj$report()
r2 <- m2$tmb_obj$report()
r3 <- m3$tmb_obj$report()

## -----------------------------------------------------------------------------
sd3$Estimate[row.names(sd3) == "ar1_phi"]
sd3$Estimate[row.names(sd3) == "ar1_phi"] +
  c(-2, 2) * sd3$`Std. Error`[row.names(sd3) == "ar1_phi"]

## -----------------------------------------------------------------------------
plot_map_raster <- function(dat, column = "est") {
  ggplot(dat, aes_string("X", "Y", fill = column)) +
    geom_raster() +
    facet_wrap(~year) +
    coord_fixed() +
    scale_fill_viridis_c()
}

## -----------------------------------------------------------------------------
p1 <- predict(m1, newdata = qcs_grid)
p2 <- predict(m2, newdata = qcs_grid)
p3 <- predict(m3, newdata = qcs_grid)

## -----------------------------------------------------------------------------
plot_map_raster(filter(p1, year == 2003), "zeta_s")
plot_map_raster(filter(p2, year == 2003), "zeta_s")
plot_map_raster(filter(p3, year == 2003), "zeta_s")

## -----------------------------------------------------------------------------
plot_map_raster(p1, "est")
plot_map_raster(p2, "est")
plot_map_raster(p3, "est")

## -----------------------------------------------------------------------------
plot_map_raster(p2, "est_rf") + scale_fill_gradient2()
plot_map_raster(p3, "est_rf") + scale_fill_gradient2()

## -----------------------------------------------------------------------------
plot_map_raster(filter(p1, year == 2003), "omega_s")
plot_map_raster(filter(p2, year == 2003), "omega_s")
plot_map_raster(filter(p3, year == 2003), "omega_s")


d <- readRDS("~/src/gfsynopsis/report/data-cache/pacific-cod.rds")
# d <- readRDS("~/src/gfsynopsis/report/data-cache/pacific-ocean-perch.rds")
# d <- readRDS("~/src/gfsynopsis/report/data-cache/silvergray-rockfish.rds")
# d <- readRDS("~/src/gfsynopsis/report/data-cache/redbanded-rockfish.rds")
d <- d$survey_sets
# d <- dplyr::filter(d, survey_series_id == 1)
# d <- gfplot::get_survey_sets("pacific cod", ssid = 1)


library(dplyr)
# dat <- gfplot:::tidy_survey_sets(d, "SYN WCVI", years = seq(1, 1e6))
dat <- gfplot:::tidy_survey_sets(d, "SYN QCS", years = seq(1, 1e6))
# dat <- gfplot:::tidy_survey_sets(dat, "HBLL OUT S", years = seq(1, 3000), density_column = "density_ppkm2")
# d <- d %>% select(year, longitude, latitude, depth_m, density_kgpm2)

dat <- mutate(dat, density = density*1000*1000)
dat <- filter(dat, !is.na(depth))
dat <- gfplot:::scale_survey_predictors(dat)
dat <- select(dat, -X10, -Y10)

pcod <- dat
usethis::use_data(pcod, internal = FALSE, overwrite = TRUE)

grid_locs <- gfplot:::make_prediction_grid(filter(dat, year == 2003), survey = "SYN QCS", cell_width = 2)$grid
# grid_locs <- gfplot:::make_prediction_grid(filter(dat, year == 2004), survey = "SYN WCVI", cell_width = 2)$grid
grid_locs <- rename(grid_locs, depth = akima_depth)
grid_locs$year <- NULL
qcs_grid <- grid_locs
usethis::use_data(qcs_grid, internal = FALSE, overwrite = TRUE)

# pcod_spde <- make_spde(pcod$X, pcod$Y, n_knots = 140)
# plot_spde(pcod_spde)
# m <- sdmTMB(
#  pcod, density ~ 0 + as.factor(year) + depth_scaled + depth_scaled2,
#  time = "year", spde = pcod_spde, family = tweedie(link = "log"), silent = FALSE,
#   anisotropy = TRUE
# )
# plot_anisotropy(m)
# 
# pcod$resids <- residuals(m) # randomized quantile residuals
# hist(pcod$resids)
# qqnorm(pcod$resids);abline(a = 0, b = 1)
# 
# library(ggplot2)
# ggplot(pcod, aes(X, Y, col = resids)) + scale_colour_gradient2() +
#   geom_point() + facet_wrap(~year)
# 
# 
# predictions <- predict(m, newdata = grid_locs)
# 
# plot_map <- function(dat, column) {
#   ggplot(dat, aes_string("X", "Y", fill = column)) +
#     geom_raster() +
#     facet_wrap(~year) +
#     coord_fixed()
# }
# 
# plot_map(predictions$data, "exp(est)") +
#   scale_fill_viridis_c(trans = "sqrt") +
#   ggtitle("Prediction (fixed effects + all random effects)")
# 
# plot_map(predictions$data, "exp(est_fe)") +
#   ggtitle("Prediction (fixed effects only)") +
#   scale_fill_viridis_c(trans = "sqrt")
# 
# plot_map(predictions$data, "est_re_s") +
#   ggtitle("Spatial random effects only") +
#   scale_fill_gradient2()
# 
# plot_map(predictions$data, "est_re_st") +
#   ggtitle("Spatiotemporal random effects only") +
#   scale_fill_gradient2()
# 
# ind <- get_index(predictions, bias_correct = FALSE) # not bias correcting for speed
# 
# library(ggplot2)
# scale <- 2*2/(1000*1000)
# ggplot(ind, aes(year, est*scale)) + geom_line() +
#   geom_ribbon(aes(ymin = lwr*scale, ymax = upr*scale), alpha = 0.4)
# 
# knitr::kable(ind, format = "pandoc")
# 
# 
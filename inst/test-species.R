library(dplyr)
library(sdmTMB)
d <- readRDS("~/src/gfsynopsis/report/data-cache/buffalo-sculpin.rds")
d <- readRDS("~/src/gfsynopsis/report/data-cache/canary-rockfish.rds")
d <- d$survey_sets

dat <- gfplot:::tidy_survey_sets(d, "SYN QCS", years = 2017)

dat <- mutate(dat, density = density*1000*1000)
dat <- filter(dat, !is.na(depth))
dat <- gfplot:::scale_survey_predictors(dat)
dat <- select(dat, -X10, -Y10)

grid_locs <- gfplot:::make_prediction_grid(filter(dat, year == 2017), survey = "SYN QCS", cell_width = 2)$grid
grid_locs <- rename(grid_locs, depth = akima_depth)
grid_locs$year <- NULL

spde <- make_spde(dat$X, dat$Y, n_knots = 200)
plot_spde(spde)

xmat <- poly(log(dat$depth), 3)
coefs <- attr(xmat, "coefs")
dat$dp <- xmat[,1]
dat$dp2 <- xmat[,2]
dat$dp3 <- xmat[,3]

m <- sdmTMB(
 data = dat, formula = density ~ dp + dp2 + dp3,
 time = "year", spde = spde, family = tweedie(link = "log"), silent = FALSE,
 anisotropy = FALSE
)
# plot_anisotropy(m)

x <- dat$dp
x2 <- dat$dp2
x3 <- dat$dp3
idx <- order(x)
y <- m$model$par[[1]] + x * m$model$par[[2]] + x2 * m$model$par[[3]] + x3 * m$model$par[[4]]
plot(dat$depth[idx], exp(y[idx]), type = "o", log = "x")

dat$resids <- residuals(m) # randomized quantile residuals
hist(dat$resids)
qqnorm(dat$resids);abline(a = 0, b = 1)

library(ggplot2)
ggplot(dat, aes(X, Y, col = resids)) + scale_colour_gradient2() +
  geom_point() + facet_wrap(~year)

new_poly <- poly(log(grid_locs$depth), 3, coefs = coefs)
grid_locs$dp <- new_poly[,1]
grid_locs$dp2 <- new_poly[,2]
grid_locs$dp3 <- new_poly[,3]

x <- grid_locs$dp
x2 <- grid_locs$dp2
x3 <- grid_locs$dp3
idx <- order(x)
y <- m$model$par[[1]] + x * m$model$par[[2]] + x2 * m$model$par[[3]] + x3 * m$model$par[[4]]
plot(grid_locs$depth[idx], exp(y[idx]), type = "l", log = "x")
abline(v = min(dat$depth[dat$present == 1]), col = "red")
abline(v = max(dat$depth[dat$present == 1]), col = "red")

predictions <- predict(m, newdata = grid_locs)

plot_map <- function(dat, column) {
  ggplot(dat, aes_string("X", "Y", fill = column)) +
    geom_raster() +
    facet_wrap(~year) +
    coord_fixed()
}

predictions$data <- filter(predictions$data,
  depth_scaled >= min(dat$depth_scaled[dat$present %in% c(0,1)]), depth_scaled <= max(dat$depth_scaled %in% c(0,1))
  )

plot_map(predictions$data, "exp(est)") +
  scale_fill_viridis_c(trans = "sqrt") +
  ggtitle("Prediction (fixed effects + all random effects)")

plot_map(predictions$data, "exp(est_fe)") +
  ggtitle("Prediction (fixed effects only)") +
  scale_fill_viridis_c(trans = "sqrt")

plot_map(predictions$data, "est_re_s") +
  ggtitle("Spatial random effects only") +
  scale_fill_gradient2()

plot_map(predictions$data, "est_re_st") +
  ggtitle("Spatiotemporal random effects only") +
  scale_fill_gradient2()
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

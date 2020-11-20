######################
# Estimate Centre of Gravity
######################

library(sdmTMB)
library(ggplot2)
library(dplyr)


pcod_spde <- make_spde(pcod$X, pcod$Y, n_knots = 100)

m <- sdmTMB(
  pcod, density ~ 0 + as.factor(year) + depth_scaled + depth_scaled2,
  time = "year", spde = pcod_spde, family = tweedie(link = "log"),
  silent = FALSE
)


# Predictions onto new data:
# expand for time units:

newdata <- qcs_grid # use bathymetry grid of Queen Charlotte Sound
# ggplot(qcs_grid, aes(X, Y, colour=log(depth))) + geom_point() + scale_color_viridis_c(direction= -1)

original_time <- sort(unique(m$data[[m$time]])) # select and order all unique time periods

nd <- do.call("rbind",
  replicate(length(original_time), newdata, simplify = FALSE)) # repeat all rows of 'newdata' (grid of points in bathymetry layer) once for each year
nd[[m$time]] <- rep(original_time, each = nrow(newdata)) # label rows of newdata (now 'nd') with year ('time')

p <- predict(m, newdata = nd) # use model 'm' to predict values for each point



# get_index(p, value_name = "link_total", bias_correct = FALSE)


predictions <- p$data
# A short function for plotting our predictions:
plot_map <- function(dat, column = "est") {
  ggplot(dat, aes_string("X", "Y", fill = column)) +
    geom_raster() +
    facet_wrap(~year) +
    coord_fixed()
}

plot_map(predictions, "exp(est)") +
  scale_fill_viridis_c(trans = "sqrt") +
  ggtitle("Prediction (fixed effects + all random effects)")



# Calculate centre of gravity for latitude and longitude

cog <- get_cog(p) # calculate centre of gravity for each data point
cog

ggplot(cog, aes(year, est, ymin = lwr, ymax = upr)) +
  geom_ribbon(alpha = 0.5) + geom_line() + facet_wrap(~coord, scales = "free_y") # what is free_y?


# table of COG by latitude
data.frame(Y = p$data$Y, est = exp(p$data$est), year = p$data$year) %>%
  group_by(year) %>% summarize(cog = sum(Y * est) / sum(est))






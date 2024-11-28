library(sdmTMB)
library(ggplot2)

# example with a year of all non-zeros in a delta model:
d <- pcod
d <- d[!(d$year == 2013 & d$density == 0),]
mesh <- make_mesh(d, c("X", "Y"), cutoff = 15)

# set up map and fixed values:
yrs <- unique(d$year)
.map <- seq_along(yrs)
.map[yrs == 2013] <- NA
.map <- factor(.map)
.start <- rep(0, length(yrs))
.start[yrs == 2013] <- 20

# here we only want to apply the map to the first linear predictor:
fit <- sdmTMB(
  density ~ 0 + factor(year),
  data = d,
  mesh = mesh,
  time = "year",
  spatial = "on",
  spatiotemporal = "off", # make the example fast
  control =
    sdmTMBcontrol(
      map = list(b_j = .map),
      start = list(b_j = .start)
    ),
  family = delta_gamma()
)
fit

nd <- replicate_df(qcs_grid, "year", unique(pcod$year))
pred <- predict(fit, newdata = nd, return_tmb_object = TRUE)
ind <- get_index(pred, bias_correct = TRUE)

ggplot(ind, aes(year, est, ymin = lwr, ymax = upr)) +
  geom_ribbon(fill = "grey50") +
  geom_line() +
  geom_vline(xintercept = 2013, lty = 2)

# what about Poisson-link?

# https://pbs-assess.github.io/sdmTMB/articles/poisson-link.html
# I would think you need to map the first linear predictor
# which predicts log(numbers) to a large value
# since p = 1 - exp(-n)
# so same as before
fit <- sdmTMB(
  density ~ 0 + factor(year),
  data = d,
  mesh = mesh,
  time = "year",
  spatial = "on",
  spatiotemporal = "off", # make the example fast
  control =
    sdmTMBcontrol(
      map = list(b_j = .map),
      start = list(b_j = .start)
    ),
  family = delta_gamma(type = "poisson-link")
)

fit$sd_report
# but that induces weird 2nd linear predictors since r = nw / p
# so, p is effectively fixed at 1, therefore:
# r = n * w, but n = fixed high, so weight must be low...

# so... I'm not sure what we'd want to do in this case
# thankfully no zeros is probably a pretty rare case
# in contexts where you'd use the poisson-link delta models

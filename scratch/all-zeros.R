library(sdmTMB)
library(ggplot2)

# example with a year of all zeros in a delta model:
d <- pcod
mesh <- make_mesh(d, c("X", "Y"), cutoff = 15)
d$density[d$year == 2013] <- 0

# set up map and fixed values:
yrs <- unique(d$year)
.map <- seq_along(yrs)
.map[yrs == 2013] <- NA
.map <- factor(.map)
.start <- rep(0, length(yrs))
.start[yrs == 2013] <- -20

# debatable whether we also want to adjust weights?
# .weights <- rep(1, nrow(d))
# .weights[d$year == 2013] <- 0

fit <- sdmTMB(
  density ~ 0 + factor(year),
  data = d,
  mesh = mesh,
  time = "year",
  # weights = .weights,
  spatiotemporal = "off", # make the example fast
  control =
    sdmTMBcontrol(
      map = list(b_j = .map, b_j2 = .map),
      start = list(b_j = .start, b_j2 = .start)
    ),
  family = delta_gamma()
)
fit

# the sanity warning is about the NA in the standard error; ignore
# maybe we can pass a flag to avoid triggering sanity() here?

# note the log likelihood will be affected by the fixed value
# this is only an issue if comparing models to other distributions
# where this parameter isn't fixed (e.g., the Tweedie)
# also need to think about how many 'parameters' we have
logLik(fit)

nd <- replicate_df(qcs_grid, "year", unique(pcod$year))
pred <- predict(fit, newdata = nd, return_tmb_object = TRUE)
ind <- get_index(pred, bias_correct = TRUE)
ind$est[ind$year == 2013] <- 0
ind$se[ind$year == 2013] <- NA

ggplot(ind, aes(year, est, ymin = lwr, ymax = upr)) +
  geom_ribbon(fill = "grey50") +
  geom_line()

# what about Poisson-link?

fit <- sdmTMB(
  density ~ 0 + factor(year),
  data = d,
  mesh = mesh,
  time = "year",
  # weights = .weights,
  spatiotemporal = "off", # make the example fast
  control =
    sdmTMBcontrol(
      map = list(b_j = .map, b_j2 = .map),
      start = list(b_j = .start, b_j2 = .start)
    ),
  family = delta_gamma(type = "poisson-link")
)
fit

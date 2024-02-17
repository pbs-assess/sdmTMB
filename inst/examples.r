d <- readRDS("inst/eulachon qPCR 2019 and 2021 standards clean.rds")
d$Ct[which(is.na(d$Ct))] = 0 # delta-models in sdmTMB need this

library(sdmTMB)

# Build a mesh to implement the SPDE approach:
mesh <- make_mesh(pcod_2011, c("X", "Y"), cutoff = 20)
fit <- sdmTMB(
  density ~ s(depth),
  data = pcod_2011, mesh = mesh,
  family = tweedie(link = "log")
)

# Use the same model as above, but pass in a 2nd dataframe for std curve estimation
# fit2 <- sdmTMB( # should throw error, because wrong family
#   density ~ s(depth),
#   data = pcod_2011, mesh = mesh,
#   family = tweedie(link = "log"),
#   control = sdmTMBcontrol(stdcurve = TRUE, stdcurve_df = d)
# )
#
# fit3 <- sdmTMB(
#   density ~ s(depth),
#   data = pcod_2011, mesh = mesh,
#   family = delta_gamma(),
#   control = sdmTMBcontrol(stdcurve_df = d)
# )

# load in observed data
d_obs <- readRDS("inst/eulachon qPCR 2019 and 2021 samples clean.rds")
d_obs$Ct[which(is.na(d_obs$Ct))] = 0 # delta-models in sdmTMB need this
d_obs$utm.lon.m <- d_obs$utm.lon.m / 1000
d_obs$utm.lat.m <- d_obs$utm.lat.m / 1000
mesh <- make_mesh(d_obs, c("utm.lon.m", "utm.lat.m"), cutoff = 20)

# Fitting the model with no standard curve (intercept only)
fit4 <- sdmTMB(
  Ct ~ 1,
  data = d_obs, mesh = mesh,
  family = delta_gaussian()
)

# Fit everything together -- intercept only. AIC 15213.08
fit5 <- sdmTMB(
  Ct ~ 1,
  data = d_obs, mesh = mesh,
  control = sdmTMBcontrol(stdcurve_df = d)
)
#sanity(fit5) # pass

#adding year: doesn't converge
fit6 <- sdmTMB(
  Ct ~ year,
  data = d_obs, mesh = mesh,
  family = delta_gaussian(),
  control = sdmTMBcontrol(stdcurve_df = d)
)

fit7 <- sdmTMB(
  Ct ~ year,
  data = d_obs, mesh = mesh,
  family = delta_gaussian(),
  time = "year",
  spatiotemporal = "iid",
  control = sdmTMBcontrol(stdcurve_df = d)
)


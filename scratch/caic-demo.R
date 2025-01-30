
library(sdmTMB)

# Build a mesh to implement the SPDE approach:
mesh <- make_mesh(pcod_2011, c("X", "Y"), cutoff = 20)

# Fit a Tweedie spatial random field GLMM with a smoother for depth:
fit <- sdmTMB(
  density ~ s(depth),
  data = pcod_2011, mesh = mesh,
  family = tweedie(link = "log"),
  control = sdmTMBcontrol(profile="b_j")
)

CAIC.sdmTMB(fit, what="CAIC")
CAIC.sdmTMB(fit, what="EDF")

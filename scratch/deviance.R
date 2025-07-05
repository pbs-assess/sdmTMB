# library(sdmTMB)
# mesh <- make_mesh(pcod_2011, c("X", "Y"), cutoff = 20)
# fit <- sdmTMB(
#   density ~ depth_scaled,
#   data = pcod_2011, mesh = mesh,
#   family = tweedie(link = "log")
# )
# fit
# r <- fit$tmb_obj$report()


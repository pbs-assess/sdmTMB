set.seed(123)
# a1 is a fake predictor:
predictor_dat <- data.frame(
  X = runif(300), Y = runif(300),
  a1 = rnorm(300), year = rep(1, each = 50)
)
mesh <- make_mesh(predictor_dat, xy_cols = c("X", "Y"), cutoff = 0.1)

sim_dat <- sdmTMB_sim2(
  formula = ~ 1 + a1,
  data = predictor_dat,
  time = "year",
  mesh = mesh,
  family = gaussian(link = "identity"),
  range = 0.5,
  phi = 0.05,
  sigma_O = 0.25,
  seed = 3542,
  B = c(0.2, -0.4) # B0 = intercept, B1 = a1 slope
)

# m <- sdmTMB(observed ~ a1, data = sim_dat, mesh = mesh, experimental = list(sim_re = TRUE))
# x <- as.list(m$sd_report, "Estimate")$omega_s
# s <- m$tmb_obj$simulate(m$tmb_obj$env$last.par.best)
# # plot(x, s$omega_s)
# plot(s$y_i, sim_dat$observed)

m <- sdmTMB(observed ~ a1, data = sim_dat, mesh = mesh, do_fit = FALSE,
  experimental = list(sim_re = TRUE))

tmb_params <- m$tmb_params
tmb_params$b_j <- c(0, 10)
tmb_params$ln_phi <- log(20)
tmb_params$ln_tau_O <- 1

obj <- TMB::MakeADFun(m$tmb_data, parameters = tmb_params, random = m$tmb_random,
  map = m$tmb_map, DLL = "sdmTMB")
s <- obj$simulate()

plot(s$omega_s)

# plot(s$y_i)
plot(m$tmb_data$X_ij[,2], s$y_i)


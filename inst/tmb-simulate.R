set.seed(1)
data <- data.frame(
  X = runif(300), Y = runif(300),
  a1 = rnorm(300), year = rep(1, each = 50)
)
mesh <- make_mesh(predictor_dat, xy_cols = c("X", "Y"), cutoff = 0.1)
formula = ~ 1 + a1
response <- sdmTMB:::get_response(formula)
if (length(response) == 0L) {
  formula <- as.formula(paste("sdmTMB_response_", paste(as.character(formula), collapse = "")))
  data[["sdmTMB_response_"]] <- 0.1 # fake! does nothing but lets sdmTMB parse the formula
}
# get tmb_data structure; parsed model matrices etc.:
fit <- sdmTMB(
  formula = formula, data = data, mesh = mesh, time = NULL,
  family = gaussian(link = "identity"), time_varying = NULL, do_fit = FALSE,
  experimental = list(sim_re = TRUE)
)
p <- predict(fit)
obj <- fit$tmb_obj
params <- fit$tmb_params

range <- 0.5
sigma_O <- 0.2
kappa <- sqrt(8)/range
tau <- 1/(sqrt(4 * pi) * kappa * sigma_O)

params$ln_tau_O <- log(tau)
params$ln_kappa <- rep(log(kappa), 2)
params$b_j <- c(0.2, -0.4)
params$ln_phi <- log(0.1)

newobj <- TMB::MakeADFun(data = fit$tmb_data, map = fit$tmb_map,
  random = fit$tmb_random, parameters = params, DLL = "sdmTMB")

s <- newobj$simulate()

s2 <- sdmTMB_sim2(
  formula = ~ 1 + a1,
  data = predictor_dat,
  time = "year",
  mesh = mesh,
  family = gaussian(link = "identity"),
  range = 0.5,
  sigma_E = 0,
  phi = 0.1,
  sigma_O = 0.2,
  seed = 3542,
  B = c(0.2, -0.4) # B0 = intercept, B1 = a1 slope
)

data$sim_obs <- s$y_i
data$sim_omega_s <- s$omega_s_A

library(ggplot2)
ggplot(data, aes(X, Y, colour = sim_omega_s)) + geom_point() +
  scale_color_gradient2()

ggplot(data, aes(X, Y, colour = sim_obs)) + geom_point() +
  scale_color_gradient2()

ggplot(s2, aes(X, Y, colour = observed)) + geom_point() +
  scale_color_gradient2()

ggplot(s2, aes(X, Y, colour = omega_s)) + geom_point() +
  scale_color_gradient2()


m1 <- sdmTMB(observed ~ a1, data = s2, mesh = mesh)
m2 <- sdmTMB(sim_obs ~ a1, data = data, mesh = mesh)

m1
m2

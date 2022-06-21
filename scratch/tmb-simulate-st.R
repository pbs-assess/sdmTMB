set.seed(1)
data <- data.frame(
  X = runif(3000), Y = runif(3000),
  a1 = rnorm(3000), year = rep(1:10, each = 300)
)
mesh <- make_mesh(data, xy_cols = c("X", "Y"), cutoff = 0.1)
formula = ~ 1 + a1
response <- sdmTMB:::get_response(formula)
if (length(response) == 0L) {
  formula <- as.formula(paste("sdmTMB_response_", paste(as.character(formula), collapse = "")))
  data[["sdmTMB_response_"]] <- 0.1 # fake! does nothing but lets sdmTMB parse the formula
}
# get tmb_data structure; parsed model matrices etc.:
fit <- sdmTMB(
  formula = formula, data = data, mesh = mesh, time = "year",
  family = gaussian(link = "identity"), time_varying = NULL, do_fit = FALSE,
  experimental = list(sim_re = TRUE)
)
obj <- fit$tmb_obj
params <- fit$tmb_params

range <- 0.5
sigma_O <- 0.2
sigma_E <- 0.1
kappa <- sqrt(8)/range
tau_O <- 1/(sqrt(4 * pi) * kappa * sigma_O)
tau_E <- 1/(sqrt(4 * pi) * kappa * sigma_E)

params$ln_tau_O <- log(tau_O)
params$ln_tau_E <- log(tau_E)
params$ln_kappa <- rep(log(kappa), 2)
params$b_j <- c(0.2, -0.4)
params$ln_phi <- log(0.1)

newobj <- TMB::MakeADFun(data = fit$tmb_data, map = fit$tmb_map,
  random = fit$tmb_random, parameters = params, DLL = "sdmTMB")

set.seed(3928)
s <- newobj$simulate()

# check against a pure R version:
s2 <- sdmTMB_sim2(
  formula = ~ 1 + a1,
  data = data,
  time = "year",
  mesh = mesh,
  family = gaussian(link = "identity"),
  range = 0.5,
  sigma_E = 0.2,
  phi = 0.1,
  sigma_O = 0.2,
  seed = 3542,
  B = c(0.2, -0.4) # B0 = intercept, B1 = a1 slope
)

data$sim_obs <- s$y_i
data$sim_omega_s <- s$omega_s_A
data$sim_epsilon_st <- s$epsilon_st_A_vec

library(ggplot2)
ggplot(data, aes(X, Y, colour = sim_omega_s)) + geom_point() +
  scale_color_gradient2()

ggplot(s2, aes(X, Y, colour = omega_s)) + geom_point() +
  scale_color_gradient2()

ggplot(data, aes(X, Y, colour = sim_epsilon_st)) + geom_point() +
  scale_color_gradient2() + facet_wrap(~year)

ggplot(s2, aes(X, Y, colour = epsilon_st)) + geom_point() +
  scale_color_gradient2() + facet_wrap(~year)

ggplot(data, aes(X, Y, colour = sim_obs)) + geom_point() +
  scale_color_gradient2() + facet_wrap(~year)

ggplot(s2, aes(X, Y, colour = observed)) + geom_point() +
  scale_color_gradient2() + facet_wrap(~year)

m1 <- sdmTMB(observed ~ a1, data = s2, mesh = mesh, time = "year")
m2 <- sdmTMB(sim_obs ~ a1, data = data, mesh = mesh, time = "year")

m1
m2

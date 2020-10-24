library(dplyr)
library(testthat)
library(ggplot2)
library(sdmTMB)
library(future)
library(assertthat)

plan(multisession, workers = floor(availableCores() / 2))
options(future.rng.onMisuse = "ignore")

SEED <- 42
set.seed(SEED)
x <- runif(200, -1, 1)
y <- runif(200, -1, 1)
N <- length(x)
time_steps <- 10

sigma_E <- 0.2
phi <- 0.1
.range <- 0.5
b_epsilon_logit = 0.2
b_epsilon = 2*plogis(b_epsilon_logit) - 1
sigma_O = 0
rho = 0
X <- model.matrix(~ -1+as.factor(x1), data.frame(x1 = sort(rep(1:time_steps,N))))
betas = runif(10, 0, 2)
sigma_E_vec = exp(log(sigma_E) + b_epsilon * seq(0, time_steps-1))

true <- tribble(
  ~variable, ~true_value,
  "sigma_E.est", sigma_E,
  "range.est", .range,
  "b_est", b_epsilon_logit
)

loc <- data.frame(x = x, y = y)
spde <- make_mesh(loc, xy_cols = c("x", "y"), cutoff = 0.1)
plot(spde)

#
# local.stcov <- function(coords, kappa.s, variance = 1, nu = 1) {
#   s <- as.matrix(dist(coords))
#   scorr <- exp((1 - nu) * log(2) + nu * log(s * kappa.s) -
#       lgamma(nu)) * besselK(s * kappa.s, nu)
#   diag(scorr) <- 1
#   return(variance * scorr)
# }

out <- furrr::future_map(seq_len(80), function(i) {

  dat <- sdmTMB_sim(
    x = x, y = y, mesh = spde, time_steps = time_steps, rho = rho,
    betas = betas, X = X, phi = phi, range = .range,
    sigma_O = sigma_O, sigma_E = sigma_E_vec,
    seed = SEED * i, family = gaussian()
  )

  mesh <- make_mesh(dat, xy_cols = c("x", "y"), cutoff = 0.1)

  # fit log-linear model with year fixed effects
  m <- sdmTMB(
    data = dat, formula = observed ~ -1 + as.factor(time),
    time = "time", spde = mesh, reml = TRUE, include_spatial = FALSE,
    epsilon_model = "loglinear"
  )

  est <- tidy(m, conf.int = TRUE)
  est_ran <- tidy(m, "ran_pars", conf.int = TRUE)

  data.frame(
    b_est = m$sd_report$par.fixed["b_epsilon_logit"],
    b_est.lwr = m$sd_report$par.fixed["b_epsilon_logit"] - 1.96 * sqrt(m$sd_report$cov.fixed["b_epsilon_logit","b_epsilon_logit"]),
    b_est.upr = m$sd_report$par.fixed["b_epsilon_logit"] + 1.96 * sqrt(m$sd_report$cov.fixed["b_epsilon_logit","b_epsilon_logit"]),
    sigma_E.est = est_ran[est_ran$term == "sigma_E",][est_ran$term == "sigma_E", "estimate"][1],
    sigma_E.lwr = est_ran[est_ran$term == "sigma_E",][est_ran$term == "sigma_E", "estimate"][1] -
       1.96 * est_ran[est_ran$term == "sigma_E",][est_ran$term == "sigma_E", "std.error"][1],
    sigma_E.upr = est_ran[est_ran$term == "sigma_E",][est_ran$term == "sigma_E", "estimate"][1] +
       1.96 * est_ran[est_ran$term == "sigma_E",][est_ran$term == "sigma_E", "std.error"][1],
    range.est = est_ran[est_ran$term == "range","estimate"],
    range.lwr = est_ran[est_ran$term == "range","estimate"] - 1.96 * est_ran[est_ran$term == "range","std.error"],
    range.upr = est_ran[est_ran$term == "range","estimate"] + 1.96 * est_ran[est_ran$term == "range","std.error"]
    )
})

est <- bind_rows(out)
reshape2::melt(est) %>%
  right_join(true) %>%
  ggplot(aes(value)) +
  facet_wrap(vars(variable), scales = "free") +
  geom_histogram(bins = 10) +
  geom_vline(aes(xintercept = true_value), colour = "red")

coverage <- est %>%
  mutate(covered = b_est.lwr < b_est & b_est.upr > b_est) %>%
  pull(covered) %>%
  mean()
coverage


plan(sequential)

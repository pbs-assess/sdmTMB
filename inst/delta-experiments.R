library(sdmTMB)
library(ggplot2)

pcod_spde <- make_mesh(pcod, c("X", "Y"), cutoff = 15)
pcod_pos <- subset(pcod, density > 0)
pcod_spde_pos <- make_mesh(pcod_pos, c("X", "Y"), mesh = pcod_spde$mesh)

fit_tw <- sdmTMB(density ~ 1,
  data = pcod, mesh = pcod_spde,
  time = "year", family = tweedie()
)
fit_tw$sd_report

p1 <- predict(fit_tw, newdata = pcod)
p2 <- predict(fit_tw, newdata = NULL)
head(p1)
head(p2)
expect_equal(p1$est, p2$est)
expect_equal(p1$omega_s, p2$omega_s)
expect_equal(p1$epsilon_st, p2$epsilon_st)
expect_equal(p1$est_rf, p2$est_rf)
expect_equal(p1$est_non_rf, p2$est_non_rf)

p <- predict(fit_tw, newdata = qcs_grid, return_tmb_object = TRUE)
system.time(ind_tw1 <- get_index(p, bias_correct = FALSE))
system.time(ind_tw2 <- get_index(p, bias_correct = TRUE))
system.time(ind_tw3 <- get_index(p, bias_correct = TRUE))
head(ind_tw3)

ggplot(ind_tw2, aes(year, est, ymin = lwr, ymax = upr)) +
  geom_ribbon(alpha = 0.5) +
  geom_line() +
  geom_ribbon(data = ind_tw3, alpha = 0.5, fill = "red") +
  geom_line(data = ind_tw3, col = "red")

fit_dg <- sdmTMB(density ~ 1,
  data = pcod, mesh = pcod_spde,
  time = "year", family = delta_gamma(link1 = "logit", link2 = "log")
)
fit_dg$sd_report
p <- predict(fit_dg, newdata = qcs_grid)
head(p)
p <- predict(fit_dg, newdata = qcs_grid, type = "response")
head(p)

p <- predict(fit_dg, newdata = qcs_grid, return_tmb_object = TRUE)
ind_dg <- get_index(p, bias_correct = FALSE)
head(ind_dg)

ind_dgb <- get_index(p, bias_correct = TRUE)
head(ind_dgb)

# check
fit_bin <- sdmTMB(present ~ 1,
  data = pcod, mesh = pcod_spde,
  time = "year", family = binomial()
)
fit_bin$sd_report

# check
fit_gamma <- sdmTMB(density ~ 1,
  data = pcod_pos, mesh = pcod_spde_pos,
  time = "year", family = Gamma(link = "log")
)
fit_gamma$sd_report

pbin <- predict(fit_bin, newdata = qcs_grid, nsim = 500)
ppos <- predict(fit_gamma, newdata = qcs_grid, nsim = 500)
ind_sims <- get_index_sims(log(plogis(pbin) * exp(ppos)))

ggplot(ind_sims, aes(year, est, ymin = lwr, ymax = upr)) +
  geom_ribbon(alpha = 0.5) +
  geom_line() +
  geom_ribbon(data = ind_dg, alpha = 0.5, fill = "red") +
  geom_line(data = ind_dg, col = "red")

fit_dln <- sdmTMB(density ~ 1,
  data = pcod, mesh = pcod_spde,
  time = "year", family = delta_lognormal()
)
fit_dln$sd_report

# check
fit_ln <- sdmTMB(density ~ 1,
  data = pcod_pos, mesh = pcod_spde_pos,
  time = "year", family = lognormal(link = "log")
)
fit_ln$sd_report

fit_plg <- sdmTMB(density ~ 1,
  data = pcod, mesh = pcod_spde,
  time = "year", family = delta_poisson_link_gamma()
)
fit_plg$sd_report

p <- predict(fit_plg, newdata = qcs_grid, type = "response")
head(p)

p <- predict(fit_plg, newdata = pcod, type = "response")
head(p)
plot(log(p$est + 1), log(p$density + 1));abline(0, 1)

expect_error(p <- predict(fit_plg, newdata = NULL, type = "response"))

pcod$count <- round(pcod$density)

fit_dtnb2 <- sdmTMB(count ~ 1,
  data = pcod, mesh = pcod_spde,
  time = "year", family = delta_truncated_nbinom2()
)
fit_dtnb2$sd_report

pcod_spde <- make_mesh(pcod, c("X", "Y"), cutoff = 15)
pcod$offset <- rnorm(nrow(pcod))
fit1 <- sdmTMB(density ~ 1 + offset,
  data = pcod, mesh = pcod_spde,
  time = "year", family = tweedie()
)
fit2 <- sdmTMB(density ~ 1,
  offset = pcod$offset,
  data = pcod, mesh = pcod_spde,
  time = "year", family = tweedie()
)
fit3 <- sdmTMB(density ~ 1,
  data = pcod, mesh = pcod_spde,
  time = "year", family = tweedie()
)
expect_error(fit4 <- sdmTMB(density ~ 1 + offset,
  data = pcod, mesh = pcod_spde,
  time = "year", family = delta_gamma()
))
fit1$sd_report
fit2$sd_report
fit3$sd_report

# smooths
pcod_spde <- make_mesh(pcod, c("X", "Y"), cutoff = 15)
fit1 <- sdmTMB(present ~ s(depth),
  data = pcod, mesh = pcod_spde,
  family = binomial()
)
plot_smooth(fit1)
fit1$sd_report

fit2 <- sdmTMB(density ~ s(depth),
  data = pcod, mesh = pcod_spde,
  family = delta_gamma()
)
fit2$sd_report

fit_dg <- sdmTMB(density ~ depth_scaled,
  data = pcod, mesh = pcod_spde,
  time = "year", family = delta_gamma(link1 = "logit", link2 = "log")
)
fit_tw <- sdmTMB(density ~ depth_scaled,
  data = pcod, mesh = pcod_spde,
  time = "year", family = tweedie()
)

print_model_info(fit_dg)
print_model_info(fit_tw)


print_main_effects(fit_tw)
print_main_effects(fit_dg)
print_main_effects(fit_dg, 2)

fit_dg <- sdmTMB(density ~ s(depth_scaled),
  data = pcod, mesh = pcod_spde,
  family = delta_gamma(link1 = "logit", link2 = "log")
)
fit_tw <- sdmTMB(density ~ s(depth_scaled),
  data = pcod, mesh = pcod_spde,
  family = tweedie()
)

print_smooth_effects(fit_tw)
print_smooth_effects(fit_dg, 1)
print_smooth_effects(fit_dg, 2)

pcod$fyear <- as.factor(pcod$year)
fit_dg <- sdmTMB(density ~ 1 + (1 | fyear),
  data = pcod, mesh = pcod_spde,
  family = delta_gamma(link1 = "logit", link2 = "log")
)
fit_tw <- sdmTMB(density ~ 1 + (1 | fyear),
  data = pcod, mesh = pcod_spde,
  family = tweedie()
)

print_iid_re(fit_dg, 1)
print_iid_re(fit_dg, 2)
print_iid_re(fit_tw)

fit_dg <- sdmTMB(density ~ 1,
  time_varying = ~ 0 + depth_scaled,
  data = pcod, mesh = pcod_spde,
  time = "year", family = delta_gamma(link1 = "logit", link2 = "log")
)
fit_tw <- sdmTMB(density ~ 1 + (1 | fyear),
  time_varying = ~ 0 + depth_scaled,
  data = pcod, mesh = pcod_spde,
  time = "year", family = tweedie()
)

print_time_varying(fit_tw)
print_time_varying(fit_dg)
print_time_varying(fit_dg, 2)

print_range(fit_tw)
print_range(fit_dg)
print_range(fit_dg, 2)

fit_dg <- sdmTMB(density ~ 1,
  data = pcod, mesh = pcod_spde,
  time = "year", share_range = FALSE,
  family = delta_gamma(link1 = "logit", link2 = "log")
)
fit_tw <- sdmTMB(density ~ 1,
  data = pcod, mesh = pcod_spde,
  time = "year", share_range = FALSE,
  family = tweedie()
)

print_range(fit_tw)
print_range(fit_dg)
print_range(fit_dg, 2)

fit_dg <- sdmTMB(density ~ s(log(depth), k = 5),
  data = pcod, mesh = pcod_spde,
  time = "year",
  family = delta_gamma(link1 = "logit", link2 = "log")
)
fit_tw <- sdmTMB(density ~ s(log(depth), k = 5),
  data = pcod, mesh = pcod_spde,
  time = "year",
  family = tweedie()
)

print_other_parameters(fit_tw)
print_other_parameters(fit_dg)
print_other_parameters(fit_dg, 2)

print_one_model(fit_tw)
print_one_model(fit_dg, 1)
print_one_model(fit_dg, 2)

print(fit_dg)
print(fit_tw)

# --------------

#
#
# # here we will fix the random field parameters at their approximate
# # MLEs (maximum likelihood estimates) from a previous fit
# # to improve speed of convergence:
# m_tmb <- sdmTMB(present ~ 0 + as.factor(year),
#   data = pcod, mesh = pcod_spde, family = binomial(),
#   time = "year")
# m_tmb
#
#
#
# p <- m_tmb$model$par
# ln_kappa <- p[names(p) == "ln_kappa"][[1]]
# ln_tau_O <- p[names(p) == "ln_tau_O"][[1]]
# ln_tau_E <- p[names(p) == "ln_tau_E"][[1]]
#
# m_tmb2 <- sdmTMB(present ~ 0 + as.factor(year),
#   data = pcod, mesh = pcod_spde, family = binomial(),
#   time = "year",
#   control = sdmTMBcontrol(
#     start = list(ln_kappa = rep(ln_kappa, 2),
#     ln_tau_E = ln_tau_E, ln_tau_O = ln_tau_O),
#     map = list(ln_kappa = rep(factor(NA), 2),
#       ln_tau_E = factor(NA), ln_tau_O = factor(NA)))
# )
#
# system.time({
#   m_tmb3 <- sdmTMB(present ~ 0 + as.factor(year),
#     data = pcod, mesh = pcod_spde, family = binomial(),
#     time = "year",
#     priors = sdmTMBpriors(b = normal(rep(0, 9), rep(2, 9))),
#     control = sdmTMBcontrol(
#       start = list(ln_kappa = rep(ln_kappa, 2)),
#       map = list(ln_kappa = rep(factor(NA), 2)))
#   )
# })
#
# library(tmbstan)
# options(mc.cores = parallel::detectCores())
#
# system.time({
#   m_stan3 <- tmbstan(m_tmb3$tmb_obj, iter = 600L, warmup = 200L,
#     chains = 4L, thin = 1L
#   )
# })
# print(m_stan3, pars = c("b_j", "omega_s[1]", "epsilon_st[1]", "lp__", "ln_tau_E", "ln_tau_O"))
# print(m_stan3, pars = c("b_j", "ln_tau_E", "ln_tau_O"))
#
#
# library(tmbstan)
# options(mc.cores = parallel::detectCores())
# m_stan2 <- tmbstan(m_tmb2$tmb_obj, iter = 400L, warmup = 150L,
#   chains = 4L, thin = 1L, init = "last.par.best"
# )
# print(m_stan2, pars = c("b_j", "omega_s[1]", "epsilon_st[1]", "lp__"))
#
# pred <- predict(m_tmb2, tmbstan_model = m_stan2, newdata = qcs_grid)
#
# ind <- get_index_sims(pred)
#
# ggplot(ind, aes(year, est)) + geom_line() +
#   geom_ribbon(aes(ymin = lwr, ymax = upr), alpha = 0.4)
#
# e <- extract(m_stan3)
#
#

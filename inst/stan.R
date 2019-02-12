library(sdmTMB)
# d <- subset(pcod, year >= 2009) # subset for example speed
d <- pcod
pcod_spde <- make_spde(d$X, d$Y, n_knots = 60)
plot_spde(pcod_spde)

library(ggplot2)
ggplot(d, aes(X, Y, size = density)) + geom_point() +
  facet_wrap(~year)

# d <- subset(d, density > 0)
d$log_density <- log(d$density)
nrow(d)

ggplot(d, aes(X, Y, size = density)) + geom_point() +
  facet_wrap(~year)

m <- sdmTMB(
  d, density ~ 1 + as.factor(year) + depth_scaled + depth_scaled2,
  time = "year", spde = pcod_spde, family = tweedie(link = "log"),
  silent = FALSE, enable_priors = TRUE)

library(tmbstan)
options(mc.cores = parallel::detectCores())
m$tmb_obj$retape()

chains <- 4L
inits <- lapply(seq_len(chains), function(x) m$model$par)

fit <- tmbstan(m$tmb_obj, chains = chains, iter = 200,
  lower = c(-1, rep(-5, length(m$model$par)-1)),
  upper = c(10, rep(5, length(m$model$par)-1)),
  init = inits, laplace = TRUE)

e <- rstan::extract(fit)
names(e)
m$model$par

traceplot(fit, pars = c("b_j", "ln_tau_O", "ln_kappa", "ln_phi"))
plot(density(exp(e$ln_tau_O)))
plot(density(exp(e$ln_kappa)))
plot(density(e$b_j[,10]))
plot(density(e$b_j[,11]))

x_pred <- seq(min(d$depth_scaled), max(d$depth_scaled), length.out = 300)
b0 <- e$b_j[,1] # pick a year
b10 <- e$b_j[,10]
b11 <- e$b_j[,11]

pred <- purrr::map_df(sample(seq_len(nrow(e$b_j)), 200L), function(i) {
  tibble::tibble(
    i = i,
    y_hat = exp(b0[i] + b10[i] * x_pred + b11[i] * x_pred^2),
    x = -exp((x_pred * d$depth_sd[1] + d$depth_mean[1])))
})

tmb_pred <- data.frame(x = -exp((x_pred * d$depth_sd[1] + d$depth_mean[1])),
  y_hat = exp(m$model$par[[1]] + m$model$par[[10]] * x_pred + m$model$par[[11]] * x_pred^2))

ggplot(pred, aes(x, y_hat, group = i)) + geom_line(alpha = 0.05) +
  coord_cartesian(xlim = c(-400, 0), expand = FALSE,
    ylim = c(0, max(pred$y_hat) * 1.05)) +
  xlab("Depth (m)") +
  ylab("Predicted density in some units") +
  geom_line(data = tmb_pred, aes(x, y_hat), inherit.aes = FALSE, col = "red") +
  gfplot::theme_pbs()

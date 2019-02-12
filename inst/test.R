library(ggplot2)
library(doParallel)
library(dplyr)
library(sdmTMB)

cores <- parallel::detectCores()
registerDoParallel(cores = cores)

# -----------------------------------------------------------------------------
# Check the parameter estimates:
out <- plyr::ldply(seq_len(cores * 1), function(i) {
  dat <- sim(
    time_steps = 6,
    plot = FALSE, seed = i*1,
    sigma_O = 0.2,
    sigma_E = 0.3,
    phi = 0.4, kappa = 1.3)
  spde <- make_spde(x = dat$x, y = dat$y, n_knots = 100)
  # plot_spde(spde)

  m <- sdmTMB(data = dat, formula = z ~ 1, time = "time",
    family = gaussian(link = "identity"), spde = spde)
  r <- m$tmb_obj$report()
  tibble::tibble(
    sigma_O = r$sigma_O,
    sigma_E = r$sigma_E,
    kappa = exp(r$ln_kappa),
    phi = exp(r$ln_phi)
  )

  # group_by(dat, time) %>%
    # summarise(total = sum(exp(z))) %>%
    # mutate(i = i)
  # predict(m)$data %>% mutate(i = i)

  # predict(m, newdata = expand.grid(x = seq(0, 10, length.out = )))
  # s <- summary(TMB::sdreport(m$tmb_obj))
  # head(summary(s))
}, .parallel = TRUE)

out

tidyr::gather(out) %>%
  dplyr::left_join(tibble::tibble(
    key = c("sigma_O", "sigma_E", "kappa", "phi"),
    true_value = c(0.2, 0.3, 1.3, 0.4)
  )) %>%
  ggplot(aes(key, value)) + geom_point() +
  geom_point(aes(y = true_value), colour = "red")

# -----------------------------------------------------------------------------
# Check the spatial predictions:

out <- plyr::ldply(seq_len(cores * 1), function(i) {
  dat <- sim(
    time_steps = 6,
    plot = FALSE, seed = i*1,
    sigma_O = 0.2,
    sigma_E = 0.3,
    phi = 0.4, kappa = 1.3)
  spde <- make_spde(x = dat$x, y = dat$y, n_knots = 100)
  m <- sdmTMB(data = dat, formula = z ~ 1, time = "time",
    family = gaussian(link = "identity"), spde = spde)
  r <- m$tmb_obj$report()
  predict(m)$data %>% mutate(i = i)
}, .parallel = TRUE)

ggplot(filter(out, i %in% 1:3), aes(real_z, est)) +
  geom_point() +
  facet_grid(time~i) +
  coord_equal() +
  geom_abline(colour = "red")

ggplot(filter(out, i %in% 1:3), aes(est_re_s, omega_s)) +
  geom_point() +
  facet_grid(i~time) +
  coord_equal() +
  geom_abline(colour = "red")

ggplot(filter(out, i %in% 1:3), aes(est_re_st, epsilon_st)) +
  geom_point() +
  facet_grid(i~time) +
  coord_equal() +
  geom_abline(colour = "red")

ggplot(filter(out, i %in% 1:3, time == 1),
  aes(x, y, colour = est_re_s)) + geom_point() +
  facet_grid(~i) +
  coord_equal() +
  scale_color_gradient2()

ggplot(filter(out, i %in% 1:3, time == 1),
  aes(x, y, colour = omega_s)) + geom_point() +
  facet_grid(~i) +
  coord_equal() +
  scale_color_gradient2()

ggplot(filter(out, i %in% 1:3, time == 1:3),
  aes(x, y, colour = est_re_st)) + geom_point() +
  facet_grid(time~i) +
  coord_equal() +
  scale_color_gradient2()

ggplot(filter(out, i %in% 1:3, time %in% 1:3),
  aes(x, y, colour = epsilon_st)) + geom_point() +
  facet_grid(time~i) +
  coord_equal() +
  scale_color_gradient2()

# ----------------------------------------------------------------
# Check the projections to new data

# FIXME: CRASH!


out <- plyr::llply(seq_len(1), function(i) {
  dat <- sim(
    time_steps = 6,
    plot = FALSE, seed = i*1,
    sigma_O = 0.2,
    sigma_E = 0.3,
    phi = 0.4, kappa = 1.3)
  dat$X <- dat$x
  dat$Y <- dat$y
  spde <- make_spde(x = dat$X, y = dat$Y, n_knots = 200)
  m <- sdmTMB(data = dat, formula = z ~ 1, time = "time",
    family = gaussian(link = "identity"), spde = spde)
  r <- m$tmb_obj$report()
  nd <- expand.grid(X = seq(0, 10, 0.05),
    Y = seq(0, 10, 0.05))
  nd$time <- 1
  list(
    data = dat,
    pred = predict(m, newdata = nd)$data %>% mutate(i = i))
}, .parallel = FALSE)

g <- ggplot(out[[1]]$pred, aes(X, Y, fill = est)) +
  geom_raster() +
  facet_grid(~time) +
  coord_equal() +
  scale_fill_viridis_c()

g1 <- ggplot(out[[1]]$data, aes(X, Y, colour = real_z)) +
  geom_point(cex = 1.2) +
  facet_grid(~time) +
  coord_equal() +
  scale_colour_viridis_c()

gridExtra::grid.arrange(g, g1, nrow = 2)

# -----------------------------------------------------------------
# SE

dat <- sim(
  x = runif(100),
  y = runif(100),
  time_steps = 1,
  plot = FALSE, seed = sample.int(1000, 1),
  sigma_O = 0.8,
  sigma_E = 0.3,
  phi = 0.2, kappa = 1.2)
dat$X <- dat$x
dat$Y <- dat$y
spde <- make_spde(x = dat$X, y = dat$Y, n_knots = 95)
plot_spde(spde)
m <- sdmTMB(data = dat, formula = z ~ 1, time = "time",
  family = gaussian(link = "identity"), spde = spde)
# nd <- expand.grid(X = seq(0, 1, 0.2),
#   Y = seq(0, 1, 0.2))
# nd$time <- 1
pred <- predict(m, newdata = dat, se_fit = TRUE)$data

pred$i <- seq_len(nrow(pred))
library(ggplot2)
ggplot(pred, aes(i, est, ymin = est - 2 * est_se, ymax = est + 2 * est_se)) +
  geom_linerange() +
  geom_point(aes(y = real_z)) +
  coord_flip()

mean((pred$real_z < pred$est - 1.96 * pred$est_se) |
    (pred$real_z > pred$est + 1.96 * pred$est_se))

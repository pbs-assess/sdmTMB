library(ggplot2)
library(dplyr)
devtools::load_all()

set.seed(123)

s <- purrr::map_dfr(1:3, function(i) {
  predictor_dat <- data.frame(
    X = runif(300), Y = runif(300),
    year = rep(1:6, each = 50)
  )
  mesh <- make_mesh(predictor_dat, xy_cols = c("X", "Y"), cutoff = 0.1)
  sdmTMB_simulate(
    formula = ~1,
    data = predictor_dat,
    time = "year",
    mesh = mesh,
    family = gaussian(),
    range = 0.5,
    sigma_E = 0.2,
    phi = 0.1,
    sigma_O = 0,
    spatiotemporal = "rw",
    seed = i,
    B = 0
  )
}, .id = "category")
row.names(s) <- NULL
head(s)

ggplot(s, aes(X, Y, colour = epsilon_st)) +
  geom_point() +
  facet_grid(category ~ year) +
  scale_colour_viridis_c()

mesh <- make_mesh(s, c("X", "Y"), cutoff = 0.1)
fit <- sdmTMB(
  observed ~ 0 + as.factor(category),
  data = s,
  time = "year",
  group = "category",
  spatial = "off",
  spatiotemporal = "rw",
  silent = FALSE,
  mesh = mesh,
  do_fit = T
)

# fit$tmb_obj$par
# fit$tmb_params$upsilon_stc
# dim(fit$tmb_params$upsilon_stc)
#
# fit$tmb_map$ln_kappa
# fit$tmb_map$epsilon_st
# fit$tmb_map$epsilon_re
# fit$tmb_map$upsilon_stc
# fit$tmb_map$ln_tau_E
# fit$tmb_map$ln_tau_O
# fit$tmb_data$n_c
# fit$tmb_data$mvrw_cat_i
# fit$tmb_data$year_i
# fit$tmb_data$rw_fields
# fit$tmb_data$ar1_fields

nd <- expand.grid(X = seq(0, 1, length.out = 50), Y = seq(0, 1, length.out = 50), category = unique(s$category), year = unique(s$year))
p1 <- predict(fit, newdata = NULL)
p2 <- predict(fit, newdata = s)
plot(p1$est, p2$est)

plot(p2$est, s$mu)


p_temp <- predict(fit, newdata = nd, return_tmb_report = TRUE)

p <- predict(fit, newdata = nd)
p$upsilon_stc <- p_temp$proj_upsilon_st_A_vec

est <- as.list(fit$sd_report, "Estimate", report = TRUE)
se <- as.list(fit$sd_report, "Std. Error", report = TRUE)
est$log_sigma_U
est$sigma_U
se$log_sigma_U

ggplot(p, aes(X, Y, colour = upsilon_stc)) +
  geom_point() +
  facet_grid(category ~ year) +
  scale_colour_gradient2()

ggplot(s, aes(X, Y, colour = mu)) +
  geom_point() +
  facet_grid(category ~ year) +
  scale_colour_viridis_c()

ggplot(p, aes(X, Y, fill = est)) +
  geom_raster() +
  facet_grid(category ~ year) +
  scale_fill_viridis_c()

d <- pcod
pcod_spde <- make_spde(d$X, d$Y, n_knots = 75)
plot_spde(pcod_spde)

# Tweedie:
m_rw1 <- sdmTMB(d, density ~ 0 + as.factor(year),
  time = "year", time_varying = ~ 0 + depth_scaled + depth_scaled2,
  spde = pcod_spde, family = tweedie(link = "log"),
  silent = FALSE, ar1_fields = FALSE, include_spatial = TRUE)

m_rw1$tmb_obj$report()$ln_tau_V
r <- m_rw1$tmb_obj$report()
r$b_j
r$b_rw_t


get_y_hat <- function(b0, b1, b2, year) {
  x_pred <- seq(min(d$depth_scaled), max(d$depth_scaled), length.out = 300)
  data.frame(
    x = -exp((x_pred * d$depth_sd[1] + d$depth_mean[1])),
    y_hat = exp(b0 + b1 * x_pred + b2 * x_pred^2),
    year = year)
}

n_t <- nrow(r$b_rw_t)
yrs <- sort(unique(d$year))
pred <- purrr::map_df(seq_len(n_t), function(.t) {
  get_y_hat(r$b_j[.t], r$b_rw_t[.t,1], r$b_rw_t[.t,2], yrs[.t])
})

library(ggplot2)
ggplot(pred, aes(x = x, ymax = y_hat + year * 5, ymin = year * 5,
  group = year, fill = year)) +
  geom_ribbon(lwd = 0.3, alpha = 0.5, col = "grey50") +
  xlab("Depth (m)") +
  scale_color_viridis_c(option = "C") +
  scale_fill_viridis_c(option = "C") +
  ylab("Predicted density in some units") +
  gfplot::theme_pbs() +
  coord_flip(xlim = c(0, -350))

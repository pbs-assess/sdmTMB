library(dplyr)
library(testthat)
library(ggplot2)
library(sdmTMB)
library(future)
plan(multisession, workers = floor(availableCores() / 2))
options(future.rng.onMisuse = "ignore")

SEED <- 42
set.seed(SEED)
x <- runif(200, -1, 1)
y <- runif(200, -1, 1)
N <- length(x)
time_steps <- 9
betas <- c(0.5, 0.7)
sigma_E <- 0.3
phi <- 0.2
rho <- 0.5
.range <- 0.8
X <- model.matrix(~ x1, data.frame(x1 = rnorm(N * time_steps)))

true <- tribble(
  ~variable, ~true_value,
  "sigma_E.est", sigma_E,
  "b0.est", betas[1],
  "b1.est", betas[2],
  "range.est", .range,
  "phi.est", phi,
  "rho.est", rho
)

loc <- data.frame(x = x, y = y)
spde <- make_mesh(loc, xy_cols = c("x", "y"), cutoff = 0.1)
plot(spde)

out <- furrr::future_map(seq_len(8*10), function(i) {
  s <- sdmTMB_sim(
    x = x, y = y, mesh = spde, X = X, sigma_V = c(0, 0),
    betas = betas, time_steps = time_steps, ar1_phi = rho,
    phi = phi, range = .range, sigma_O = 0, sigma_E = sigma_E,
    seed = SEED * i, family = gaussian()
  )
  # ggplot(s, aes(x, y, colour = mu)) +
  #   geom_point() +
  #   scale_colour_gradient2() +
  #   facet_wrap(vars(time))
  mesh <- make_mesh(s, xy_cols = c("x", "y"), cutoff = 0.1)
  m <- sdmTMB(
    data = s, formula = observed ~ x1,
    time = "time", spde = mesh, reml = TRUE,
    ar1_fields = TRUE, include_spatial = FALSE
  )
  est <- tidy(m, conf.int = TRUE)
  est_ran <- tidy(m, "ran_pars", conf.int = TRUE)
  data.frame(
    b0.est = est[est$term == "(Intercept)", "estimate"],
    b0.lwr = est[est$term == "(Intercept)", "conf.low"],
    b0.upr = est[est$term == "(Intercept)", "conf.high"],
    b1.est = est[est$term == "x1", "estimate"],
    b1.lwr = est[est$term == "x1", "conf.low"],
    b1.upr = est[est$term == "x1", "conf.high"],
    rho.est = est_ran[est_ran$term == "rho", "estimate"],
    rho.lwr = est_ran[est_ran$term == "rho", "conf.low"],
    rho.upr = est_ran[est_ran$term == "rho", "conf.high"],
    sigma_E.est = est_ran[est_ran$term == "sigma_E", "estimate"],
    sigma_E.lwr = est_ran[est_ran$term == "sigma_E", "conf.low"],
    sigma_E.upr = est_ran[est_ran$term == "sigma_E", "conf.high"],
    phi.est = est_ran[est_ran$term == "phi", "estimate"],
    phi.lwr = est_ran[est_ran$term == "phi", "conf.low"],
    phi.upr = est_ran[est_ran$term == "phi", "conf.high"],
    range.est = est_ran[est_ran$term == "range", "estimate"],
    range.lwr = est_ran[est_ran$term == "range", "conf.low"],
    range.upr = est_ran[est_ran$term == "range", "conf.high"]
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
  mutate(covered = sigma_E.lwr < sigma_E & sigma_E.upr > sigma_E) %>%
  pull(covered) %>%
  mean()
coverage

coverage <- est %>%
  mutate(covered = b0.lwr < betas[1] & b0.upr > betas[1]) %>%
  pull(covered) %>%
  mean()
coverage

worm_plot <- function(dat, lwr, est, upr, true, title) {
  median_est <- pull(dat, {{est}}) %>% median()
  dat %>%
    arrange({{est}}) %>%
    mutate(sim = seq_len(n())) %>%
    mutate(covered = {{lwr}} < true & {{upr}} > true) %>%
    ggplot(aes({{est}}, sim, xmin = {{lwr}}, xmax = {{upr}}, shape = covered)) +
    geom_vline(xintercept = true, col = "red", lty = 2) +
    geom_vline(xintercept = median_est, col = "grey50", lty = 1) +
    geom_point(size = 1.5) +
    geom_linerange(size = 0.3) +
    scale_shape_manual(values = c("TRUE" = 19, "FALSE" = 21)) +
    guides(shape = FALSE) +
    ggtitle(title) +
    labs(x = "Parameter value") +
    theme(axis.title.y = element_blank(), axis.text.y = element_blank(),
      plot.title = element_text(hjust = 0.5))
}
theme_set(ggsidekick::theme_sleek())
g1 <- worm_plot(est, rho.lwr, rho.est, rho.upr, rho, expression(rho))
g2 <- worm_plot(est, b0.lwr, b0.est, b0.upr, betas[1], expression(beta[0]))
g3 <- worm_plot(est, b1.lwr, b1.est, b1.upr, betas[2], expression(beta[1]))
g3 <- worm_plot(est, phi.lwr, phi.est, phi.upr, phi, expression(phi))
g4 <- worm_plot(est, sigma_E.lwr, sigma_E.est, sigma_E.upr, sigma_E, expression(sigma[Epsilon]))
g5 <- worm_plot(est, range.lwr, range.est, range.upr, .range, expression(range))
cowplot::plot_grid(g1, g2, g3, g4, g5, ncol = 3)

coverage <- est %>%
  mutate(covered = b1.lwr < betas[2] & b1.upr > betas[2]) %>%
  pull(covered) %>%
  mean()
coverage

coverage <- est %>%
  mutate(covered = rho.lwr < rho.true & rho.upr > rho.true) %>%
  pull(covered) %>%
  mean()
coverage

# expect_equal(median(est$b0.est), betas[1], tol = 0.01)
# expect_equal(median(est$phi), phi, tol = 0.01)
# expect_equal(median(est$sigma_E), sigma_E, tol = 0.01)
# expect_equal(median(est$range), .range, tol = 0.01)
# expect_gt(coverage, 0.92)
# expect_lt(coverage, 0.98)

plan(sequential)

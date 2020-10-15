library(dplyr)
library(testthat)
library(ggplot2)
library(sdmTMB)
library(future)
plan(multisession, workers = floor(availableCores() / 2) - 1)
options(future.rng.onMisuse = "ignore")

SEED <- 42
set.seed(SEED)
x <- stats::runif(200, -1, 1)
y <- stats::runif(200, -1, 1)
time_steps <- 9
betas <- 0.5
sigma_E <- 0.3
phi <- 0.2
ar1_phi <- 0.5
.range <- 1.5

true <- tibble(
  variable = c("sigma_E", "b.est", "range", "phi", "rho"),
  true_value = c(sigma_E, betas, .range, phi, ar1_phi)
)

loc <- data.frame(x = x, y = y)
spde <- make_mesh(loc, xy_cols = c("x", "y"), n_knots = 130)
plot(spde)

out <- furrr::future_map(seq(1, 7 * 10), function(i) {
# out <- map(seq(1, 1), function(i) {
  s <- sdmTMB_sim(
    x = x, y = y, mesh = spde,
    betas = betas, time_steps = time_steps, ar1_phi = ar1_phi,
    phi = phi, range = .range, sigma_O = 0, sigma_E = sigma_E,
    seed = SEED * i * 3
  )
  # ggplot(s, aes(x, y, colour = mu)) +
  #   geom_point() +
  #   scale_colour_gradient2() +
  #   facet_wrap(vars(time))
  mesh <- make_mesh(s, xy_cols = c("x", "y"), n_knots = 130)
  # plot(mesh)
  m <- tryCatch(
    {
      sdmTMB(
        data = s, formula = observed ~ 0 + cov1,
        time = "time", spde = spde, reml = TRUE,
        ar1_fields = TRUE, include_spatial = FALSE
      )
    },
    error = function(e) NA
  )
  if (identical(m, NA)) {
    return(NA)
  }
  est <- tidy(m, conf.int = TRUE)
  est_ran <- tidy(m, "ran_pars", conf.int = TRUE)
  p <- as.list(m$model$par)
  r <- m$tmb_obj$report()
  data.frame(
    sigma_E = r$sigma_E,
    rho = r$rho,
    b.est = est[est$term == "cov1", "estimate"],
    b.lwr = est[est$term == "cov1", "conf.low"],
    b.upr = est[est$term == "cov1", "conf.high"],
    rho.lwr = est_ran[est_ran$term == "ar1_phi", "conf.low"],
    rho.upr = est_ran[est_ran$term == "ar1_phi", "conf.high"],
    range = r$range,
    phi = exp(p$ln_phi)
  )
})

d <- bind_rows(out)
reshape2::melt(d) %>%
  left_join(true) %>%
  ggplot(aes(value)) +
  facet_wrap(vars(variable), scales = "free") +
  geom_histogram(bins = 20) +
  geom_vline(aes(xintercept = true_value), colour = "red")

coverage <- d %>%
  mutate(covered = b.lwr < betas, b.upr > betas) %>%
  pull(covered) %>%
  mean()
coverage

coverage <- d %>%
  mutate(covered = rho.lwr < ar1_phi, rho.upr > ar1_phi) %>%
  pull(covered) %>%
  mean()
coverage

expect_equal(median(d$b.est), betas, tol = 0.01)
expect_equal(median(d$phi), phi, tol = 0.01)
expect_equal(median(d$sigma_E), sigma_E, tol = 0.01)
expect_equal(median(d$range), .range, tol = 0.01)
# expect_gt(coverage, 0.92)
# expect_lt(coverage, 0.98)

plan(sequential)

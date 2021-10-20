library(dplyr)
library(testthat)
library(ggplot2)
library(sdmTMB)
library(future)
plan(multisession, workers = floor(availableCores() / 2))
options(future.rng.onMisuse = "ignore")

SEED <- 1
set.seed(SEED)
x <- runif(500, -1, 1)
y <- runif(500, -1, 1)
N <- length(x)
betas <- c(0.5, 0.8)
sigma_O <- 0.3
phi <- 0.2
.range <- 0.8
X <- model.matrix(~ x1, data.frame(x1 = rnorm(N)))

true <- tribble(
  ~variable, ~true_value,
  "sigma_O", sigma_O,
  "b0.est", betas[1],
  "b1.est", betas[2],
  "range", .range,
  "phi", phi,
)

loc <- data.frame(x = x, y = y)
spde <- make_mesh(loc, xy_cols = c("x", "y"), n_knots = 100)
plot(spde)

out <- furrr::future_map(seq_len(100), function(i) {
  s <- sdmTMB_sim(
    x = x, y = y, mesh = spde,
    betas = betas, time_steps = 1,
    phi = phi, range = .range, sigma_O = sigma_O, sigma_E = 0,
    seed = SEED * i, X = X,
    family = gaussian(),
    # family = binomial(),
    # family = tweedie(),
    # family = nbinom2(),
    # family = lognormal(),
    # family = student(),
    # family = Beta(),
    # family = poisson(),
    # family = Gamma(link = "log"),
    tweedie_p = 1.5
  )
  m <- sdmTMB(
    observed ~ x1, data = s, spde = spde,
    family = gaussian()
    # family = binomial()
    # family = tweedie()
    # family = nbinom2()
    # family = lognormal()
    # family = student()
    # family = Beta() # FIXME!
    # family = poisson()
    # family =  Gamma(link = "log")
  )

  # e <- as.list(m$sd_report, "Estimate")
  # tweedie_p <- plogis(e$thetaf) + 1
  est <- tidy(m, conf.int = TRUE)
  p <- as.list(m$model$par)
  r <- m$tmb_obj$report()
  data.frame(
    sigma_O = r$sigma_O,
    b0.est = est[est$term == "(Intercept)", "estimate"],
    b0.lwr = est[est$term == "(Intercept)", "conf.low"],
    b0.upr = est[est$term == "(Intercept)", "conf.high"],
    b1.est = est[est$term == "x1", "estimate"],
    b1.lwr = est[est$term == "x1", "conf.low"],
    b1.upr = est[est$term == "x1", "conf.high"],
    range = r$range, phi = exp(p$ln_phi)
  )
})

est <- bind_rows(out)
reshape2::melt(est) %>%
  right_join(true) %>%
  ggplot(aes(value)) +
  facet_wrap(vars(variable), scales = "free") +
  # geom_density() +
  geom_histogram(bins = 15) +
  geom_histogram(aes(y = after_stat(density)), bins = 15) +
  geom_vline(aes(xintercept = true_value), colour = "red")

coverage <- est %>%
  mutate(covered = b1.lwr < betas[2], b1.upr > betas[2]) %>%
  pull(covered) %>%
  mean()
coverage

coverage <- est %>%
  mutate(covered = b0.lwr < betas[1], b0.upr > betas[1]) %>%
  pull(covered) %>%
  mean()
coverage

expect_equal(median(est$b1.est), betas[2], tol = 0.01)
expect_equal(median(est$b0.est), betas[1], tol = 0.01)
expect_equal(median(est$phi), phi, tol = 0.01)
expect_equal(median(est$sigma_O), sigma_O, tol = 0.01)
expect_equal(median(est$range), .range, tol = 0.01)
expect_gt(coverage, 0.92)
expect_lte(coverage, 1)

plan(sequential)

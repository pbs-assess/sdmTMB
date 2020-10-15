library(dplyr)
library(testthat)
library(ggplot2)
library(sdmTMB)
library(future)
plan(multisession, workers = floor(availableCores() / 2))
options(future.rng.onMisuse = "ignore")

SEED <- 1
set.seed(SEED)
x <- stats::runif(400, -1, 1)
y <- stats::runif(400, -1, 1)
betas <- 0.5
sigma_O <- 0.3
phi <- 0.2
.range <- 0.5

true <- tibble(
  variable = c("sigma_O", "sigma_E", "b.est", "range", "phi"),
  true_value = c(sigma_O, sigma_E, betas, .range, phi)
)

x <- stats::runif(200, -1, 1)
y <- stats::runif(200, -1, 1)
loc <- data.frame(x = x, y = y)
spde <- make_mesh(loc, xy_cols = c("x", "y"), n_knots = 130)

out <- furrr::future_map(seq(1, 8 * 20), function(i) {
  s <- sdmTMB_sim(
    x = x, y = y, mesh = spde,
    betas = betas, time = rep(1L, length(x)),
    phi = phi, range = .range, sigma_O = sigma_O, sigma_E = 0,
    seed = SEED * i * 4
  )
  m <- tryCatch(
    {
      sdmTMB(
        data = s, formula = observed ~ 0 + cov1,
        time = "time", spatial_only = TRUE, spde = spde, reml = FALSE
      )
    },
    error = function(e) NA
  )
  if (identical(m, NA)) {
    return(NA)
  }
  est <- tidy(m, conf.int = TRUE)
  p <- as.list(m$model$par)
  r <- m$tmb_obj$report()
  data.frame(
    sigma_O = r$sigma_O,
    b.est = est[est$term == "cov1", "estimate"],
    b.lwr = est[est$term == "cov1", "conf.low"],
    b.upr = est[est$term == "cov1", "conf.high"],
    range = r$range, phi = exp(p$ln_phi)
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

expect_equal(median(d$b.est), betas, tol = 0.001)
expect_equal(median(d$phi), phi, tol = 0.01)
expect_equal(median(d$sigma_O), sigma_O, tol = 0.01)
expect_equal(median(d$range), .range, tol = 0.01)
expect_gt(coverage, 0.92)
expect_lt(coverage, 0.98)

plan(sequential)

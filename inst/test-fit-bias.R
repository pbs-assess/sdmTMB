# context("Is there any bias in basic fits?")

library(dplyr)
library(ggplot2)
library(sdmTMB)
library(future)
plan(multisession, workers = floor(availableCores()/2))
options(future.rng.onMisuse = "ignore")

SEED <- 1
set.seed(SEED)
x <- stats::runif(100, -1, 1)
y <- stats::runif(100, -1, 1)

# test_that("sdmTMB model fit with a covariate beta", {
initial_betas <- 0.5
kappa <- 6 # decay of spatial correlation (smaller = slower decay)
sigma_O <- 0.3 # SD of spatial process
sigma_E <- 0.3 # SD of spatial process
phi <- 0.05 # observation error

set.seed(1)
out <- furrr::future_map(seq(1, 8*20), function(i) {
# out <- purrr::map(seq(1, 8*2), function(i) {

  # set.seed(SEED * i)
  x <- stats::runif(100, -1, 1)
  y <- stats::runif(100, -1, 1)
  s <- sim(
    x = x, y = y,
    initial_betas = initial_betas, time_steps = 1L,
    phi = phi, kappa = kappa, sigma_O = sigma_O,
    seed = SEED * i
  )
  spde <- make_mesh(s, xy_cols = c("x", "y"), n_knots = 90)
  # plot(spde)
  m <- tryCatch({sdmTMB(data = s, formula = observed ~ 0 + cov1,
    time = "time", spatial_only = TRUE, spde = spde)}, error = function(e) NA)
  if (identical(m, NA)) return(NA)
  est <- tidy(m, conf.int = TRUE)
  p <- as.list(m$model$par)
  r <- m$tmb_obj$report()
  data.frame(
    sigma_O = r$sigma_O,
    # sigma_E = r$sigma_E,
    b.est = p$b_j,
    b.lwr = est[est$term == "cov1", "conf.low"],
    b.upr = est[est$term == "cov1", "conf.high"],
    kappa = exp(p$ln_kappa), phi = exp(p$ln_phi))
})

d <- bind_rows(out)
true <- tibble(variable = c("sigma_O", "sigma_E", "b.est", "kappa", "phi"),
  true_value = c(sigma_O, sigma_E, initial_betas, kappa, phi))

reshape2::melt(d) %>% left_join(true) %>% ggplot(aes(value)) +
  facet_wrap(vars(variable), scales = "free") +
  geom_histogram(bins = 20) +
  geom_vline(aes(xintercept = true_value), colour = "red")

nrow(d)
mutate(d, covered = b.lwr < 0.5, b.upr > 0.5) %>%
  pull(covered) %>% mean()

plan(sequential)



# expect_equal(m$model$convergence, 0L)
# expect_equal((p$b_j - initial_betas)^2, 0, tol = 0.05)
# expect_equal((exp(p$ln_phi) - phi)^2, 0, tol = 0.05)
# expect_equal((r$sigma_O - sigma_O)^2, 0, tol = 0.05)
# expect_equal((r$sigma_E - sigma_E)^2, 0, tol = 0.05)
# expect_equal(exp(p$ln_kappa), kappa, tol = 1.1)
# p <- predict(m)
# r <- residuals(m)
# expect_equal(mean((p$est - s$observed)^2), 0, tol = 0.02)
# })

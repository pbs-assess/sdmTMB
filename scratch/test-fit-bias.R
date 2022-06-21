# context("Is there any bias in basic fits?")

library(dplyr)
library(ggplot2)
library(sdmTMB)
library(future)
plan(multisession, workers = floor(availableCores()/2))
options(future.rng.onMisuse = "ignore")

SEED <- 1
set.seed(SEED)
x <- stats::runif(400, -1, 1)
y <- stats::runif(400, -1, 1)

# test_that("sdmTMB model fit with a covariate beta", {
initial_betas <- 0.5
kappa <- 6 # decay of spatial correlation (smaller = slower decay)
sigma_O <- 0.3 # SD of spatial process
sigma_E <- 0.3 # SD of spatial process
phi <- 0.2 # observation error
.range <- 0.5

loc <- data.frame(x = x, y = y)

spde <- make_mesh(loc, xy_cols = c("x", "y"), n_knots = 200)
plot(spde)

s <- sdmTMB:::rspde2(cbind(x, y),
 range = sqrt(8) / 3,
 sigma = sigma_O, n = 2, mesh = spde$mesh,
 return.attributes = TRUE, seed = 2
)
loc$z <- s[,1]
ggplot(loc, aes(x, y, colour = z)) + geom_point() + scale_color_gradient2()

set.seed(1234)
out <- furrr::future_map(seq(1, 8*50), function(i) {
# out <- purrr::map(seq(1, 14), function(i) {
#
#   library(INLA)
#   source("inst/INLA-helpers.R")

  x <- stats::runif(200, -1, 1)
  y <- stats::runif(200, -1, 1)
  s <- sim(
    x = x, y = y,
    initial_betas = initial_betas, time_steps = 1L,
    phi = phi, range = .range, sigma_O = sigma_O,
    seed = SEED * i * 4
  )
  spde <- make_mesh(s, xy_cols = c("x", "y"), n_knots = 150)
  # plot(spde)
  m <- tryCatch({sdmTMB(data = s, formula = observed ~ 0 + cov1,
    time = "time", spatiotemporal = "off", spde = spde, reml = TRUE)}, error = function(e) NA)
  if (identical(m, NA)) return(NA)
  est <- tidy(m, conf.int = TRUE)
  p <- as.list(m$model$par)
  r <- m$tmb_obj$report()
  data.frame(
    sigma_O = r$sigma_O,
    b.est = est[est$term == "cov1", "estimate"],
    b.lwr = est[est$term == "cov1", "conf.low"],
    b.upr = est[est$term == "cov1", "conf.high"],
    range = r$range, phi = exp(p$ln_phi))
})

d <- bind_rows(out)
true <- tibble(variable = c("sigma_O", "sigma_E", "b.est", "range", "phi"),
  true_value = c(sigma_O, sigma_E, initial_betas, .range, phi))

reshape2::melt(d) %>% left_join(true) %>% ggplot(aes(value)) +
  facet_wrap(vars(variable), scales = "free") +
  geom_histogram(bins = 20) +
  geom_vline(aes(xintercept = true_value), colour = "red")

nrow(d)
mutate(d, covered = b.lwr < initial_betas, b.upr > initial_betas) %>%
  pull(covered) %>% mean()

plan(sequential)


mean(d$b.est)
median(d$b.est)
mean(d$sigma_O)
median(d$sigma_O)
median(d$phi)
mean(d$phi)


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

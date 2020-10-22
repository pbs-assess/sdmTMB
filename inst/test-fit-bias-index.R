library(dplyr)
library(testthat)
library(ggplot2)
library(sdmTMB)
library(future)
theme_set(ggsidekick::theme_sleek())
plan(multisession, workers = floor(availableCores() / 2))
options(future.rng.onMisuse = "ignore")

SEED <- 1
set.seed(SEED)
x <- seq(-1, 1, length.out = 25)
y <- seq(-1, 1, length.out = 25)
loc <- expand.grid(x = x, y = y)
x <- loc$x
y <- loc$y
N <- length(x)
time_steps <- 20
sigma_O <- 0.6
sigma_E <- 0.3
phi <- 1.2
thetaf <- 1.7 # Tweedie p
.range <- 1

loc <- data.frame(x = x, y = y)
spde <- make_mesh(loc, xy_cols = c("x", "y"), cutoff = 0.2)
plot(spde)

sim_test_index <- function(i) {
  betas <- rnorm(time_steps, 1, 1)
  dat <- data.frame(year = as.factor(rep(seq_len(time_steps), each = N)))
  s <- sdmTMB_sim(
    x = x, y = y, mesh = spde, X = model.matrix(~ 0 + year, dat),
    betas = betas, time_steps = time_steps,
    phi = phi, range = .range, sigma_O = sigma_O, sigma_E = sigma_E,
    seed = SEED * i, family = tweedie(), thetaf = thetaf
  )
  s <- mutate(s, year = dat$year)
  s_sampled <- s %>% group_by(year) %>% slice_sample(n = 25L)

  if (FALSE) {
    ggplot(s, aes(x, y, fill = mu)) +
      geom_raster() +
      geom_point(aes(size = observed), s_sampled, pch = 21, fill = NA) +
      scale_fill_viridis_c() +
      scale_size_area() +
      facet_wrap(vars(time)) +
      coord_cartesian(expand = FALSE)
  }

  mesh <- make_mesh(s_sampled, xy_cols = c("x", "y"), cutoff = 0.2)
  m <- tryCatch({sdmTMB(
    data = s_sampled,
    formula = observed ~ 0 + as.factor(year),
    time = "time",
    spde = mesh,
    reml = TRUE,
    family = tweedie(),
  )}, error = function(e) return(NULL))

  if (max(m$gradients) > 0.001) {
    m <- tryCatch({run_extra_optimization(m,
      nlminb_loops = 0, newton_steps = 1)},
      error = function(e) return(NULL))
  }

  grid <- select(s, x, y, time, year)
  p <- predict(m, newdata = grid, return_tmb_object = TRUE)
  index_bc <- get_index(p, bias_correct = TRUE)

  true_index <- group_by(s, year) %>%
    summarise(true = sum(mu), .groups = "drop") %>%
    mutate(year = as.integer(year))

  if (FALSE) {
    ggplot(index_bc, aes(time, est, ymin = lwr, ymax = upr)) +
      geom_ribbon(fill = "grey70") +
      geom_line() +
      geom_line(aes(year, true), true_index, colour = "red",
        inherit.aes = FALSE)
  }

  index_bc <- bind_cols(select(true_index, true), index_bc)
  bind_cols(data.frame(iter = i), index_bc)
}

est <- furrr::future_map_dfr(seq_len(8*3), sim_test_index)
saveRDS(est, "inst/sim-test-index.rds")
est <- readRDS("inst/sim-test-index.rds")

est %>% filter(max_gradient > 0.001)

ggplot(est, aes(time, est, ymin = lwr, ymax = upr)) +
  geom_ribbon(fill = "grey70") +
  geom_line() +
  geom_line(aes(y = true), colour = "red") +
  facet_wrap(vars(iter), scales = "free_y") #+
  # scale_y_log10()

coverage <- est %>%
  mutate(covered = lwr < true & upr > true) %>%
  summarise(coverage = mean(covered))
coverage

error <- est %>%
  mutate(error = (true - est)/true) %>%
  summarise(error = median(error))
error

plan(sequential)

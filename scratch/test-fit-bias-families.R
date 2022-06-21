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
x <- runif(800, -1, 1)
y <- runif(800, -1, 1)
N <- length(x)
time_steps <- 1
betas <- c(0.5, 0.7)
sigma_O <- 0.4
phi <- 0.5
.range <- 1
X <- model.matrix(~ x1, data.frame(x1 = rnorm(N * time_steps)))

true <- tribble(
  ~variable, ~true_value,
  "b0.est", betas[1],
  "b1.est", betas[2],
  "phi.est", phi,
  "range.est", .range,
  "sigma_O.est", sigma_O
)

loc <- data.frame(x = x, y = y)
spde <- make_mesh(loc, xy_cols = c("x", "y"), cutoff = 0.1)
plot(spde)

families <- list(
  Beta(),
  Gamma(link = "log"),
  # Gamma(link = "inverse"),
  binomial(),
  gaussian(),
  lognormal(),
  nbinom2(),
  poisson(),
  student(),
  tweedie()
)

est <- purrr::map_dfr(families, function(.fam) {
  cat(.fam$family, "\n")
  furrr::future_map_dfr(seq_len(60), function(i) {
   if (.fam$family == "Beta") {
     phi <- phi * 5 # small phi's can cause ~0.0 and ~1.0
   }
    if (.fam$link == "inverse") {
      betas[1] <- 5 # needs to keep all mu(i) > 0
    }
    s <- sdmTMB_sim(
      x = x, y = y, mesh = spde, X = X, sigma_V = c(0, 0),
      betas = betas, time_steps = time_steps,
      phi = phi, range = .range, sigma_O = sigma_O, sigma_E = 0,
      seed = SEED * i, family = .fam
    )
    # ggplot(s, aes(x, y, colour = mu)) +
    #   geom_point() +
    #   scale_colour_gradient2() +
    #   facet_wrap(vars(time))
    if (.fam$family == "Beta") { # small phi's can cause ~0.0 and ~1.0
      s$observed[s$observed > 0.9999] <- 0.9999
      s$observed[s$observed < 0.0001] <- 0.0001
    }
    m <- tryCatch({sdmTMB(
      data = s, formula = observed ~ x1,
      time = "time", spde = spde, reml = TRUE,
      spatiotemporal = "off", family = .fam
    )}, error = function(e) return(NULL))

    if (max(m$gradients) > 0.001) {
      m <- tryCatch({run_extra_optimization(m,
        nlminb_loops = 0, newton_steps = 1)},
        error = function(e) return(NULL))
    }

    est <- tidy(m, conf.int = TRUE)
    est_ran <- tidy(m, "ran_pars", conf.int = TRUE)
    est <- rbind(est, est_ran)

    data.frame(
      iter = i,
      family = .fam$family,
      link = .fam$link,
      max_gradient = max(m$gradients),
      b0.est = est[est$term == "(Intercept)", "estimate"],
      b0.lwr = est[est$term == "(Intercept)", "conf.low"],
      b0.upr = est[est$term == "(Intercept)", "conf.high"],
      b1.est = est[est$term == "x1", "estimate"],
      b1.lwr = est[est$term == "x1", "conf.low"],
      b1.upr = est[est$term == "x1", "conf.high"],
      phi.est = est_ran[est_ran$term == "phi", "estimate"],
      phi.lwr = est_ran[est_ran$term == "phi", "conf.low"],
      phi.upr = est_ran[est_ran$term == "phi", "conf.high"],
      range.est = est_ran[est_ran$term == "range", "estimate"],
      range.lwr = est_ran[est_ran$term == "range", "conf.low"],
      range.upr = est_ran[est_ran$term == "range", "conf.high"],
      sigma_O.est = est_ran[est_ran$term == "sigma_O", "estimate"],
      sigma_O.lwr = est_ran[est_ran$term == "sigma_O", "conf.low"],
      sigma_O.upr = est_ran[est_ran$term == "sigma_O", "conf.high"],
      stringsAsFactors = FALSE
    )
  })
})

saveRDS(est, "inst/sim-test-family.rds")
est <- readRDS("inst/sim-test-family.rds")
est <- est %>% mutate(
  phi.lwr = if_else(family == "Beta", phi.lwr / 5, phi.lwr),
  phi.est = if_else(family == "Beta", phi.est / 5, phi.est),
  phi.upr = if_else(family == "Beta", phi.upr / 5, phi.upr)
)

est %>% filter(max_gradient > 0.002)
est %>% filter(sigma_O.est > 2)

est <- est %>% filter(max_gradient < 0.002) %>%
  filter(sigma_O.est < 2)

est %>%
  reshape2::melt() %>%
  right_join(true, by = "variable") %>%
  filter(!(family == "binomial" & variable == "phi.est")) %>%
  filter(!(family == "poisson" & variable == "phi.est")) %>%
  ggplot(aes(value)) +
  facet_grid(vars(paste(family, link)), vars(variable), scales = "free") +
  geom_histogram(bins = 15) +
  geom_vline(aes(xintercept = true_value), colour = "red")

coverage <- est %>%
  filter(family != "poisson") %>%
  filter(family != "binomial") %>%
  group_by(family, link) %>%
  mutate(covered = phi.lwr < phi & phi.upr > phi) %>%
  summarise(coverage = mean(covered))
coverage

worm_plot <- function(dat, lwr, est, upr, true, title) {
  median_est <- pull(dat, {{est}}) %>% median()
  dat %>%
    mutate(sim = seq_len(n())) %>%
    mutate(covered = {{lwr}} < true & {{upr}} > true) %>%
    ggplot(aes({{est}}, sim, xmin = {{lwr}}, xmax = {{upr}}, shape = covered)) +
    geom_vline(xintercept = true, col = "red", lty = 2) +
    geom_point(size = 1.5) +
    geom_linerange(size = 0.3) +
    scale_shape_manual(values = c("TRUE" = 19, "FALSE" = 21)) +
    guides(shape = FALSE) +
    ggtitle(title) +
    labs(x = "Parameter value") +
    theme(axis.title.y = element_blank(), axis.text.y = element_blank(),
      plot.title = element_text(hjust = 0.5))
}

est %>%
  filter(family != "poisson") %>%
  filter(family != "binomial") %>%
  worm_plot(phi.lwr, phi.est, phi.upr, phi, expression(phi)) +
  coord_cartesian(xlim = c(phi-0.2, phi+0.2)) +
  facet_wrap(~paste(family, link), scales = "free")

plan(sequential)

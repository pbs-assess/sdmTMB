library(sdmTMB)

# match glmmTMB:
ctl <- sdmTMBcontrol(newton_loops = 0, multiphase = FALSE)

N <- 3e4
set.seed(1)
dat <- data.frame(
  y = fishMod::rTweedie(N, 1, 1, 1.5)
)

bench::mark({m2 <- sdmTMB(
  y ~ 1,
  data = dat,
  spatial = 'off',
  family = tweedie("log"),
  control = ctl
)}, iterations = 1)

tidy(m2)

bench::mark({m1 <- glmmTMB::glmmTMB(
  y ~ 1,
  data = dat,
  family = glmmTMB::tweedie("log")
)}, iterations = 1)


N <- 3e4
set.seed(1)
dat <- data.frame(
  y = rnorm(N, 0, 1)
)

bench::mark({m2 <- sdmTMB(
  y ~ 1,
  data = dat,
  spatial = 'off',
  control = sdmTMBcontrol(newton_loops = 1, multiphase = FALSE)
)}, iterations = 1)

bench::mark({m1 <- glmmTMB::glmmTMB(
  y ~ 1,
  data = dat
)}, iterations = 1)


##############

N <- 3e4
set.seed(1)
dat <- data.frame(
  y = rgamma(N, 0.2, 0.2/3)
)

bench::mark({m2 <- sdmTMB(
  y ~ 1,
  data = dat,
  spatial = 'off',
  family = Gamma(link = "log"),
  control = sdmTMBcontrol(newton_loops = 1, multiphase = FALSE)
)}, iterations = 1)

bench::mark({m1 <- glmmTMB::glmmTMB(
  y ~ 1,
  data = dat,
  family = Gamma(link = "log"),
)}, iterations = 1)

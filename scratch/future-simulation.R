library(ggplot2)

mesh <- make_mesh(pcod_2011, c("X", "Y"), cutoff = 20)

fit <- sdmTMB(
  density ~ 1,
  data = pcod_2011, mesh = mesh, time = "year", spatiotemporal = "rw",
  family = tweedie(link = "log")
)
fit

p <- get_pars(fit)
dim(p$epsilon_st)
eps <- p$epsilon_st
b <- tidy(fit, "ran_pars")

# sim just last year epsilon_st:
eps2 <- array(dim = dim(eps) + c(0, 0, 0))
eps2[,,1] <- cbind(eps[,-4,1], rep(NA, dim(eps)[1]))

s <- sdmTMB_simulate(
   ~ 1,
  data = pcod_2011,
  mesh = mesh,
  range = b$estimate[b$term == "range"],
  sigma_E = b$estimate[b$term == "sigma_E"],
  sigma_O = b$estimate[b$term == "sigma_O"],
  phi = b$estimate[b$term == "phi"],
  B = unname(coef(fit)),
  time = "year",
  spatiotemporal = "rw",
  family = tweedie(link = "log"),
  # seed = 92,
  fixed_re = list(epsilon_st = eps2)
)

ggplot(s, aes(X, Y, colour = epsilon_st)) + geom_point() +
  scale_color_viridis_c(limits = c(-3, 3)) +
  facet_wrap(~year)

pred <- predict(fit, newdata = pcod_2011)
ggplot(pred, aes(X, Y, colour = epsilon_st)) + geom_point() +
  scale_color_viridis_c(limits = c(-3, 3)) +
  facet_wrap(~year)

# library(dplyr)
# group_by(s, year) |>
#   summarise(.sd = sd(epsilon_st))
# group_by(pred, year) |>
#   summarise(.sd = sd(epsilon_st))

# future sim year:
p <- get_pars(fit)
eps <- p$epsilon_st
eps2 <- array(dim = dim(eps) + c(0, 1, 0))
eps2[,,1] <- cbind(eps[,,1], rep(NA, dim(eps)[1]))
b <- tidy(fit, "ran_pars")

# fake sampling locations:
samps <- sample(1:nrow(pcod_2011), size = 500)
nd <- data.frame(X = pcod_2011$X[samps], Y = pcod_2011$Y[samps], year = 2019L)
nd <- rbind(select(pcod_2011, X, Y, year), nd)

mesh2 <- make_mesh(nd, c("X", "Y"), mesh = mesh$mesh)
s <- sdmTMB_simulate(
   ~ 1,
  data = nd,
  mesh = mesh2,
  range = b$estimate[b$term == "range"],
  sigma_E = b$estimate[b$term == "sigma_E"],
  sigma_O = b$estimate[b$term == "sigma_O"],
  phi = b$estimate[b$term == "phi"],
  B = unname(coef(fit)),
  time = "year",
  spatiotemporal = "rw",
  family = tweedie(link = "log"),
  fixed_re = list(epsilon_st = eps2)
)

ggplot(s, aes(X, Y, colour = epsilon_st)) + geom_point() +
  scale_color_viridis_c(limits = c(-4, 4)) +
  facet_wrap(~year)

ggplot(pred, aes(X, Y, colour = epsilon_st)) + geom_point() +
  scale_color_viridis_c(limits = c(-4, 4)) +
  facet_wrap(~year)

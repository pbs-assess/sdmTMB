set.seed(1234)
x <- stats::runif(1000, -1, 1)
y <- stats::runif(1000, -1, 1)
loc <- data.frame(x = x, y = y)
spde <- make_mesh(loc, c("x", "y"), n_knots = 150, type = "kmeans")

s <- sdmTMB_sim(x = x, y = y,
  betas = 0, time = 1L,
  phi = 0.1, range = 2.2,
  sigma_O = 0.2, sigma_E = 0,
  seed = 93827, mesh = spde
)
library(ggplot2)
ggplot(s, aes(x, y, colour = mu)) + geom_point() +
  scale_colour_gradient2()
ggplot(s, aes(x, y, colour = observed)) + geom_point() +
  scale_colour_gradient2()

g <- rep(gl(60, 10), 999)
set.seed(134)
RE_vals <- rnorm(60, 0, 0.4)
RE_vals

h <- rep(gl(40, 10), 999)
set.seed(1283)
RE_vals2 <- rnorm(40, 0, 0.2)
RE_vals2

s$g <- g[1:nrow(s)]
s$h <- h[1:nrow(s)]
s$observed <- s$observed + RE_vals[s$g] + RE_vals2[s$h]

m1 <- sdmTMB(
  data = s,
  formula = observed ~ 1, spde = spde, silent = FALSE)
tidy(m1, "fixed", conf.int = TRUE)
tidy(m1, "ran_pars", conf.int = TRUE)

m <- sdmTMB(
  data = s,
  time = NULL,
  formula = observed ~ 1 + (1 | g) + (1 | h), spde = spde, silent = FALSE)
tidy(m, "fixed", conf.int = TRUE)
tidy(m, "ran_pars", conf.int = TRUE)
print(m)

b <- as.list(m$sd_report, "Estimate")
plot(c(RE_vals, RE_vals2), b$RE, col = m$tmb_data$ln_tau_G_index + 1, pch = 19)
abline(a = 0, b = 1)

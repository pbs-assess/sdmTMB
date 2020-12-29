set.seed(12235)
x <- stats::runif(2000, -1, 1)
y <- stats::runif(2000, -1, 1)
loc <- data.frame(x = x, y = y)
spde <- make_mesh(loc, c("x", "y"), n_knots = 150, type = "kmeans")

s <- sdmTMB_sim(x = x, y = y,
  betas = 0.5, time = 1L,
  phi = 0.01, range = 0.5,
  sigma_O = 0.2, sigma_E = 0,
  seed = 1, mesh = spde
)

g <- rep(gl(100, 10), 999)
set.seed(1)
RE_vals <- rnorm(100, 0, 0.4)
RE_vals

h <- rep(gl(95, 10), 999)
set.seed(1283)
RE_vals2 <- rnorm(95, 0, 0.2)
RE_vals2

# hist(s$observed)
s$g <- g[1:nrow(s)]
s$h <- h[1:nrow(s)]
s$observed <- s$observed + RE_vals[s$g] + RE_vals2[s$h]

# m <- sdmTMB(data = s, formula = observed ~ 1, spde = spde, silent = FALSE)
m <- sdmTMB(data = s,
  formula = observed ~ 1 + (1 | g) + (1 | h), spde = spde, silent = FALSE)

b <- as.list(m$sd_report, "Estimate")
b$RE

plot(c(RE_vals, RE_vals2), b$RE, col = m$tmb_data$ln_tau_G_index + 1, pch = 19)

exp(b$ln_tau_G)

set.seed(1)
x <- stats::runif(200, -1, 1)
y <- stats::runif(200, -1, 1)
loc <- data.frame(x = x, y = y)
spde <- make_mesh(loc, c("x", "y"), n_knots = 50, type = "kmeans")

s <- sdmTMB_sim(x = x, y = y,
  betas = 0.5, time = 1L,
  phi = 0.01, range = 0.5, sigma_O = 0.3, sigma_E = 0,
  seed = 1, mesh = spde
)

g <- gl(20, 10)
g
set.seed(123)
RE_vals <- rnorm(20, 0, 2)
RE_vals

hist(s$observed)
s$g <- g
s$observed <- s$observed + RE_vals[s$g]

# m <- sdmTMB(data = s, formula = observed ~ 1, spde = spde, silent = FALSE)
m <- sdmTMB(data = s, formula = observed ~ 1 + (1 | g), spde = spde, silent = FALSE)


library(INLA)
get_inla_mesh <- function (mesh, t.sub = 1:nrow(mesh$graph$tv)) {
    idx = cbind(mesh$graph$tv[t.sub, c(1:3, 1), drop = FALSE],
      NA)
    x = mesh$loc[t(idx), 1]
    y = mesh$loc[t(idx), 2]
    tv = mesh$graph$tv[t.sub, , drop = FALSE]
    idx = unique(as.vector(tv))
    vertices <- as.data.frame(mesh$loc[idx, , drop = FALSE])
    vertices <- vertices[,1:2]
    names(vertices) <- c("x", "y")

    list(mesh = tibble::tibble(x = x, y = y), vertices = vertices)
}


library(sdmTMB)
library(ggplot2)
library(dplyr)

set.seed(2)
N <- 100
obs <- data.frame(X = runif(N, 0, 1), Y = runif(N, 0, 1), type = "obs", year = 1L, stringsAsFactors = FALSE)

set.seed(18)
# set.seed(21)
d <- expand.grid(X = seq(0, 1, length.out = 80), Y = seq(0, 1, length.out = 80), year = 1L, type = "grid", stringsAsFactors = FALSE)
.rf_sim <- function(model, x, y) {
  out <- sdmTMB:::rf_sim(model, x, y)
  out - mean(out)
}



d <- rbind(d, obs)
orig_d <- d

sigma_O <- 0.3
sigma_Z <- 0.3
sigma_E <- 0.15
kappa <- 6.7
rf_omega <- RandomFields::RMmatern(nu = 1, var = sigma_O^2, scale = 1/kappa)
d$omega_s <- .rf_sim(model = rf_omega, d$X, d$Y)
rf_epsilon <- RandomFields::RMmatern(nu = 1, var = sigma_E^2, scale = 1/kappa)
d <- d %>% group_by(year) %>%
  group_split() %>%
  purrr::map_dfr(
    ~tibble(X = .x$X, Y = .x$Y, omega_s = .x$omega_s, year = .x$year,
      eps_st = .rf_sim(model = rf_epsilon, .x$X, .x$Y))) %>%
  ungroup()

d$value <- d$omega_s + d$eps_st
d$type <- orig_d$type

d_grid <- filter(d, type == "grid")
d_obs <- filter(d, type == "obs")

spde <- make_spde(d_obs$X, d_obs$Y, n_knots = 50)
# plot_spde(spde)

# get_inla_mesh(spde$mesh, main = NA, edge.color = "grey60", asp = 1)
# points(object$x, object$y, pch = 21, col = "#00000070")
# points(object$loc_centers, pch = 20, col = "red")

# omit <- c(6, 8, 17, 18, 43, 52, 70, 72, 77, 82, 94, 96, 117, 121, 131, 135, 139, 148)
# d_obs <- d_obs[-omit, , drop=FALSE]
#
# omit <- c(32, 107, 120, 125)
# d_obs <- d_obs[-omit, , drop=FALSE]

mesh <- get_inla_mesh(spde$mesh)

g <- ggplot(d_grid, aes(X, Y, fill = exp(value))) +
  geom_raster(interpolate = TRUE) +
  scale_fill_viridis_c(option = "B", trans = "log") +
  theme_void() +
  guides(fill = FALSE, size = FALSE) +
  geom_point(data = d_obs, mapping = aes(x = X, y = Y, size = exp(value)), pch = 21, fill = "grey90", alpha = 0.7, col = "grey90") +
  geom_point(data = d_obs, mapping = aes(x = X, y = Y, size = exp(value)), pch = 21, fill = NA, alpha = 0.9, col = "grey90") +
  scale_size(range = c(0.05, 8)) +
  coord_equal() +
  geom_path(data = mesh$mesh, mapping = aes(x, y), inherit.aes = FALSE, col = "white", alpha = 0.5, lwd = 0.35)
  # geom_point(data = mesh$vertices, mapping = aes(x, y), inherit.aes = FALSE, size = 2, colour = "grey60")
print(g)
ggsave("inst/rf.pdf", width = 7, height = 7, dpi = 300)

# plot(d_obs$X, d_obs$Y)
# identify(d_obs$X, d_obs$Y, plot=TRUE)
#
# omit <- c(6, 8, 17, 18, 43, 52, 70, 72, 77, 82, 94, 96, 117, 121, 131, 135, 139, 148)
#
# points(d_obs$X[-omit], d_obs$Y[-omit])

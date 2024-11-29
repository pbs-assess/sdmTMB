# https://en.wikipedia.org/wiki/Binary_search_algorithm

binary_search <- function(x, vec) {
  L <- 0
  R <- length(vec)
  while (L <= R) {
    m <- floor((L + R) / 2)
    cat("m =", m, "\n")
    if (vec[m] < x) {
      L <- m + 1
    } else if (vec[m] > x) {
      R <- m - 1
    } else {
      return(m)
    }
  }
}

binary_search(500, seq(1, 10000))

#' Binary search for an INLA mesh with a cutoff that matches the desired knots
#'
#' @param x A vector of x coordinates.
#' @param y A vector of y coordinates.
#' @param n_knots The desired number of knots.
#' @param min The minimum cutoff evaluated.
#' @param max The maximum cutoff evaluated.
#' @param length The number of cutoffs evaluated between `min` and `max` where the increments are log distributed.
#'
#' @return A mesh from [INLA::inla.mesh.create()].

binary_search_knots <- function(x, y,
                                n_knots,
                                min = 1e-4,
                                max = 1e4,
                                length = 1e6) {
  vec <- exp(seq(log(min), log(max), length.out = length))
  loc_xy <- cbind(x, y)
  L <- 0
  R <- length(vec)
  realized_knots <- Inf
  while (L <= R) {
    m <- floor((L + R) / 2)
    mesh <- INLA::inla.mesh.create(loc_xy, refine = TRUE, cutoff = vec[m])
    realized_knots <- mesh$n
    pretty_cutoff <- sprintf("%.2f", round(vec[m], 2))
    cat("cutoff =", pretty_cutoff, "| knots =", realized_knots)
    if (realized_knots > n_knots) {
      L <- m + 1
      cat(" |", cli::symbol$arrow_down, "\n")
    } else if (realized_knots < n_knots) {
      R <- m - 1
      cat(" |", cli::symbol$arrow_up, "\n")
    } else {
      cat(" |", cli::symbol$tick, "\n")
      return(mesh)
    }
  }
  cat("cutoff =", pretty_cutoff, "| knots =", mesh$n, "¯\\_(ツ)_/¯\n")
  mesh
}

library(sdmTMB)
x <- pcod$X
y <- pcod$Y


s <- binary_search_knots(x, y, 200)
d <- pcod
pcod_spde <- make_spde(d$X, d$Y, n_knots = 125) # only 50 knots for example speed
plot_spde(pcod_spde)

# Tweedie:
m <- sdmTMB(density ~ 0 + depth_scaled + depth_scaled2 + as.factor(year),
  data = d, time = "year", spde = pcod_spde, family = tweedie(link = "log"))
m


library(INLA)
library(inlabru)
library(ggplot2)

ggplot() + coord_equal() +
  geom_point(aes(x, y), data.frame(x = x, y = y), alpha = 0.5) +
  gg(m$spde$mesh)

out <- inla.mesh.assessment(m$spde$mesh,
  spatial.range = m$tmb_obj$report()$range,
  alpha = 2,
  dims = c(200, 200))
quantile(out$sd.dev, na.rm = TRUE)

# ggplot() + gg(out, aes(color = edge.len)) + coord_equal() +
#   scale_colour_viridis_c()

ggplot() + gg(out, aes(color = sd.dev)) + coord_equal() +
  scale_colour_viridis_c(limits = range(out$sd.dev, na.rm = TRUE)) +
  gg(m$spde$mesh)

s <- binary_search_knots(x, y, 125)
mesh <- make_spde(x, y, mesh = s)

m2 <- sdmTMB(density ~ 0 + depth_scaled + depth_scaled2 + as.factor(year),
  data = d, time = "year", spde = mesh, family = tweedie(link = "log"),
  silent = FALSE)
m2

out2<- inla.mesh.assessment(m2$spde$mesh,
  spatial.range = m2$tmb_obj$report()$range,
  alpha = 2,
  dims = c(200, 200))
ggplot() + gg(out2, aes(color = sd.dev)) + coord_equal() +
  scale_colour_viridis_c(limits = range(out2$sd.dev, na.rm = TRUE)) +
  gg(m2$spde$mesh)

quantile(out2$sd.dev, na.rm = TRUE)

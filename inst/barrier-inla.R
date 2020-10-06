# https://haakonbakkagit.github.io/btopic107.html
# try https://fromthebottomoftheheap.net/2016/03/27/soap-film-smoothers/

library(INLA)
library(fields)
library(rgeos)
library(viridisLite)

set.seed(2016)
set.inla.seed <- 2016

# dir.create("data/")
# download.file(url = "https://haakonbakkagit.github.io/data/WebSiteData-Archipelago.RData", destfile = "data-raw/WebSiteData-Archipelago.RData")

load(file = "data-raw/WebSiteData-Archipelago.RData")
df$y <- df$y.smelt

max.edge <- 0.6
bound.outer <- 4.6
mesh <- inla.mesh.2d(
  boundary = poly.water,
  loc = cbind(df$locx, df$locy),
  max.edge = c(1, 5) * max.edge,
  cutoff = 0.06,
  offset = c(max.edge, bound.outer)
)
plot(mesh, main = "Our mesh", lwd = 0.5)
points(df$locx, df$locy, col = "red")


A.i.s <- inla.spde.make.A(mesh, loc = cbind(df$locx, df$locy))

stk <- inla.stack(
  data = list(y = df$y, e = df$exposure),
  effects = list(
    s = 1:mesh$n,
    data.frame(m = 1, df[, 5:11]),
    # - m is the intercept
    iidx = 1:nrow(df)
  ),
  A = list(A.i.s, 1, 1),
  remove.unused = FALSE, tag = "est"
)

tl <- length(mesh$graph$tv[, 1])
# - the number of triangles in the mesh
posTri <- matrix(0, tl, 2)
for (t in 1:tl) {
  temp <- mesh$loc[mesh$graph$tv[t, ], ]
  posTri[t, ] <- colMeans(temp)[c(1, 2)]
}
posTri <- SpatialPoints(posTri)

# - compute the triangle positions
normal <- over(poly.water, posTri, returnList = T)
# - checking which mesh triangles are inside the normal area
normal <- unlist(normal)
barrier.triangles <- setdiff(1:tl, normal)

poly.barrier <- inla.barrier.polygon(mesh, barrier.triangles)
plot(poly.barrier)

barrier.model <- inla.barrier.pcmatern(mesh, barrier.triangles = barrier.triangles, prior.range = c(6, .5), prior.sigma = c(3, 0.01))

formula <- y ~ -1 + m + f(s, model = barrier.model) + f(iidx, model = "iid", hyper = hyper.iid)

par(mar = rep(.1, 4))
plot(df$locx, y = df$locy, pch = 20, asp = 1)
plot(poly.barrier, add = T, border = "black", col = "grey")

init <- c(0.044, 1.274, -0.596)

hyper.iid <- list(prec = list(prior = "pc.prec", param = c(3, 0.01)))

m <- inla(formula,
  data = inla.stack.data(stk),
  control.predictor = list(A = inla.stack.A(stk)),
  family = "poisson", E = e,
  control.inla = list(int.strategy = "eb"),
  control.mode = list(restart = T, theta = init)
)

local.plot.field <- function(field, mesh, xlim, ylim, ...) {
  stopifnot(length(field) == mesh$n)
  # - error when using the wrong mesh
  if (missing(xlim)) xlim <- poly.water@bbox[1, ]
  if (missing(ylim)) ylim <- poly.water@bbox[2, ]
  # - choose plotting region to be the same as the study area polygon
  proj <- inla.mesh.projector(mesh,
    xlim = xlim,
    ylim = ylim, dims = c(300, 300)
  )
  # - Can project from the mesh onto a 300x300 grid
  #   for plots
  field.proj <- inla.mesh.project(proj, field)
  # - Do the projection
  image.plot(list(x = proj$x, y = proj$y, z = field.proj),
    xlim = xlim, ylim = ylim, col = plasma(17), ...
  )
}

field <- m$summary.random$s$mean + m$summary.fixed["m", "mean"]
local.plot.field(field, mesh, zlim = c(-10, -1))
plot(poly.barrier, add = T, col = "grey")
points(df$locx, df$locy)

spde <- make_spde(df, c("locx", "locy"), cutoff = 0.02)
plot(spde)
m2 <- sdmTMB(y ~ 1, data = df, family = poisson(), spde = spde, silent = FALSE)


x <- seq(poly.water@bbox[1, 1], poly.water@bbox[1, 2], length.out = 200)
y <- seq(poly.water@bbox[2, 1], poly.water@bbox[2, 2], length.out = 200)
grid <- expand.grid(locx = x, locy = y)
p <- predict(m2, newdata = grid)

library(ggplot2)
ggplot(p, aes(locx, locy, fill = est)) + geom_raster() +
  scale_fill_viridis_c() + geom_point(data = df, aes(locx, locy), inherit.aes = FALSE)

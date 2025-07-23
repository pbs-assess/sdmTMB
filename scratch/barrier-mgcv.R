# example from ?mgcv::smooth.construct.so.smooth.spec

library(mgcv)
library(ggplot2)
library(sdmTMB)
library(sp)

# sim
fsb <- list(fs.boundary())
set.seed(55)
n <- 200
v <- runif(n) * 5 - 1
w <- runif(n) * 2 - 1
y <- fs.test(v, w, b = 1)
names(fsb[[1]]) <- c("v", "w")
ind <- inSide(fsb, x = v, y = w) ## remove outsiders
beta <- 0
covariate <- runif(n)
y <- y + beta * covariate + rnorm(n) * .1 ## add noise
extra_remove <- !(v > 2.5 & w < 0)
ind <- ind & extra_remove
y <- y[ind]
v <- v[ind]
w <- w[ind]
covariate <- covariate[ind]
n <- length(y)
dat <- data.frame(x = v, y = w, z = y, b = covariate)
xym <- cbind(fsb[[1]]$v, fsb[[1]]$w)

ggplot(data.frame(x = v, y = w, z = y), aes(x, y, colour = z, size = z)) +
  geom_point() +
  scale_colour_gradient2() +
  geom_path(data = as.data.frame(xym), aes(V1, V2), inherit.aes = FALSE)

loc_xy <- cbind(dat$x, dat$y)
bnd <- fmesher::fm_nonconvex_hull(loc_xy, convex = .6)
mesh <- fmesher::fm_mesh_2d(
  boundary = bnd,
  max.edge = c(0.25, .5),
  offset = c(0.05, 0.5),
  cutoff = c(0.05),
  min.angle = 20
)

spde <- make_mesh(dat, c("x", "y"), mesh = mesh)
plot(spde)

m <- sdmTMB(z ~ 1, data = dat, mesh = spde)
grid <- expand.grid(
  x = seq(min(dat$x) - 0.3, max(dat$x) + 0.3, length.out = 200),
  y = seq(min(dat$y) - 0.3, max(dat$y) + 0.3, length.out = 200), b = mean(covariate)
)

pred <- predict(m, newdata = grid)

ggplot(pred, aes(x, y, fill = est)) +
  geom_raster() +
  geom_point(data = dat, aes(x, y, size = z), inherit.aes = FALSE, pch = 21) +
  scale_fill_gradient2() +
  geom_path(data = as.data.frame(xym), aes(V1, V2), inherit.aes = FALSE) +
  coord_cartesian(expand = FALSE)

###############

# do the barrier mesh modifications by hand:
mesh <- spde$mesh

# the number of triangles in the mesh
tl <- length(mesh$graph$tv[, 1]) 
posTri <- matrix(0, tl, 2)
for (t in 1:tl) {
  temp <- mesh$loc[mesh$graph$tv[t, ], ]
  posTri[t, ] <- colMeans(temp)[c(1, 2)]
}
posTri <- sp::SpatialPoints(posTri)

xym <- cbind(fsb[[1]]$v, fsb[[1]]$w)

p <- Polygon(xym)
ps <- Polygons(list(p), 1)
sps <- SpatialPolygons(list(ps))

# compute the triangle positions
normal <- sp::over(sps, posTri, returnList = TRUE)

# checking which mesh triangles are inside the normal area
normal <- unlist(normal)
barrier.triangles <- setdiff(1:tl, normal)

if (FALSE) {
  plot(posTri@coords[normal, ])
  plot(posTri@coords[barrier.triangles, ])
}

barrier_spde <- INLAspacetime::mesh2fem.barrier(mesh, barrier.triangles = barrier.triangles)

bspde <- spde
bspde$spde_barrier <- barrier_spde
bspde$barrier_scaling <- c(1, 0.1)

mb <- sdmTMB(z ~ 1, data = dat, mesh = bspde)

pb <- predict(mb, newdata = grid)

m
mb

theme_set(theme_light())

cols <- RColorBrewer::brewer.pal(3, "BrBG")
g1 <- ggplot(pred, aes(x, y, fill = est)) +
  geom_raster() +
  geom_point(data = dat, aes(x, y, size = z), inherit.aes = FALSE, pch = 21) +
  scale_fill_gradient2(high = cols[1], low = cols[3], mid = cols[2]) +
  theme(legend.position = "none") +
  geom_path(data = as.data.frame(xym), aes(V1, V2), inherit.aes = FALSE) +
  coord_cartesian(expand = FALSE) +
  ggtitle("Regular mesh")

g2 <- ggplot(pb, aes(x, y, fill = est)) +
  geom_raster() +
  geom_point(data = dat, aes(x, y, size = z), inherit.aes = FALSE, pch = 21) +
  scale_fill_gradient2(high = cols[1], low = cols[3], mid = cols[2]) +
  theme(legend.position = "none") +
  geom_path(data = as.data.frame(xym), aes(V1, V2), inherit.aes = FALSE) +
  coord_cartesian(expand = FALSE) +
  ggtitle("Mesh with barrier")

cowplot::plot_grid(g1, g2, ncol = 1)

# from https://stat.ethz.ch/R-manual/R-patched/library/mgcv/html/smooth.construct.so.smooth.spec.html
require(mgcv)

##########################
## simple test function...
##########################

fsb <- list(fs.boundary())
nmax <- 100
## Simulate some fitting data, inside boundary...
set.seed(25)
n<-120
v <- runif(n)*5-1;w<-runif(n)*2-1
y <- fs.test(v,w,b=1)
names(fsb[[1]]) <- c("v","w")
ind <- inSide(fsb,x=v,y=w) ## remove outsiders
beta <- 0
covariate <- runif(n)
y <- y + beta * covariate + rnorm(n)*.1 ## add noise
y <- y[ind];v <- v[ind]; w <- w[ind]; covariate <- covariate[ind]
n <- length(y)

library(ggplot2)

dat <- data.frame(x = v, y = w, z = y, b = covariate)
library(sdmTMB)
# spde <- make_spde(dat, c("x", "y"), cutoff = 0.01)

xym <- cbind(fsb[[1]]$v, fsb[[1]]$w)
ggplot(data.frame(x = v, y = w, z = y), aes(x, y, colour = z, size = z)) + geom_point() +
  scale_colour_gradient2() +
  geom_path(data = as.data.frame(xym), aes(V1, V2), inherit.aes = FALSE)


loc_xy <- cbind(dat$x, dat$y)
bnd <- INLA::inla.nonconvex.hull(loc_xy, convex = .6)
mesh <- INLA::inla.mesh.2d(
  boundary = bnd,
  max.edge = c(0.3, .5),
  offset = -0.05,
  cutoff = c(0.08, 0.5),
  min.angle = 20
)
sp2 <- make_spde(dat, c("x", "y"), mesh = mesh)
plot(sp2)

spde <- sp2
plot(spde)
m <- sdmTMB(z ~ 1, data = dat, spde = spde, silent = FALSE)
grid <- expand.grid(x = seq(min(dat$x) - 0.3, max(dat$x) + 0.3, length.out = 200),
  y = seq(min(dat$y) - 0.3, max(dat$y) + 0.3, length.out = 200), b = mean(covariate))

pred <- predict(m, newdata = grid)

ggplot(pred, aes(x, y, fill = est)) + geom_raster() +
  geom_point(data = dat, aes(x, y, size = z), inherit.aes = FALSE, pch = 21) +
  scale_fill_gradient2()

###############

mesh <- spde$mesh

tl <- length(mesh$graph$tv[, 1]) # - the number of triangles in the mesh
posTri <- matrix(0, tl, 2)
for (t in 1:tl) {
  temp <- mesh$loc[mesh$graph$tv[t, ], ]
  posTri[t, ] <- colMeans(temp)[c(1, 2)]
}
posTri <- sp::SpatialPoints(posTri)

library(sp)
xym <- cbind(fsb[[1]]$v, fsb[[1]]$w)

p = Polygon(xym)
ps = Polygons(list(p),1)
sps = SpatialPolygons(list(ps))
# plot(sps)

# - compute the triangle positions
normal <- sp::over(sps, posTri, returnList = TRUE)

# - checking which mesh triangles are inside the normal area
normal <- unlist(normal)
barrier.triangles <- setdiff(1:tl, normal)

plot(posTri@coords[normal,])
plot(posTri@coords[barrier.triangles,])
# barrier.triangles <- as.numeric(normal)

library(INLA)
# poly.barrier <- inla.barrier.polygon(mesh, barrier.triangles)
# plot(poly.barrier)

barrier_spde <- inla.barrier.fem(mesh, barrier.triangles = barrier.triangles)

bspde <- spde
bspde$spde_barrier <- barrier_spde

mb <- sdmTMB(z ~ 1, data = dat, spde = bspde, silent = FALSE, barrier = TRUE,
  barrier_scaling = c(1, .05))

pb <- predict(mb, newdata = grid)

# ggplot(pb, aes(x, y, fill = est)) + geom_raster() +
#   geom_point(data = dat, aes(x, y, size = z), inherit.aes = FALSE, pch = 21) +
#   scale_fill_gradient2() +
#   geom_path(data = as.data.frame(xym), aes(V1, V2), inherit.aes = FALSE)
#
# ggplot(pred, aes(x, y, fill = est)) + geom_raster() +
#   geom_point(data = dat, aes(x, y, size = z), inherit.aes = FALSE, pch = 21) +
#   scale_fill_gradient2() +
#   geom_path(data = as.data.frame(xym), aes(V1, V2), inherit.aes = FALSE)
#
m
mb

# pb$diff <- scale(pred$omega_s, scale = FALSE) - scale(pb$omega_s, scale = FALSE)
# pb$diff <- pred$est - pb$est
#
# theme_set(theme_void())
# ggplot(pb, aes(x, y, fill = diff)) + geom_raster() +
#   scale_fill_gradient2() +
#   geom_point(data = dat, aes(x, y, size = z), inherit.aes = FALSE, pch = 21) +
#   scale_fill_gradient2() +
#   geom_path(data = as.data.frame(xym), aes(V1, V2), inherit.aes = FALSE)

cols <- RColorBrewer::brewer.pal(3, "BrBG")
g1 <- ggplot(pred, aes(x, y, fill = est)) + geom_raster() +
  geom_point(data = dat, aes(x, y, size = z), inherit.aes = FALSE, pch = 21) +
  scale_fill_gradient2(high = cols[1], low = cols[3], mid = cols[2]) +
  # theme(legend.position = "none") +
  geom_path(data = as.data.frame(xym), aes(V1, V2), inherit.aes = FALSE)

g2 <- ggplot(pb, aes(x, y, fill = est)) +
  geom_raster() +
  geom_point(data = dat, aes(x, y, size = z), inherit.aes = FALSE, pch = 21) +
  scale_fill_gradient2(high = cols[1], low = cols[3], mid = cols[2]) +
  # theme(legend.position = "none") +
  geom_path(data = as.data.frame(xym), aes(V1, V2), inherit.aes = FALSE)

cowplot::plot_grid(g1, g2)

# remotes::install_github("ropensci/rnaturalearthhires")
library(dplyr)
library(ggplot2)
library(sf)
library(sdmTMB)
theme_set(theme_minimal())

map_data <- rnaturalearth::ne_countries(
  scale = 10,
  returnclass = "sf", country = "canada"
)

d <- readRDS("~/src/gfsynopsis/report/data-cache/pacific-cod.rds")
d <- d$survey_sets
dat <- gfplot:::tidy_survey_sets(d,
  survey = c("SYN QCS", "SYN WCHG", "SYN HS", "SYN WCVI"),
  years = seq(1, 1e6))

st_bbox(map_data)

ggplot(map_data) + geom_sf()

bc_coast <- suppressMessages(
  st_crop(map_data,
    c(xmin = -134, ymin = 46, xmax = -120, ymax = 59)))

ggplot(bc_coast) + geom_sf()

crs_utm9 <- 3156

bc_coast <- st_transform(bc_coast, crs_utm9)
ggplot(bc_coast) + geom_sf()


surv <- dat %>% select(lon, lat, density) %>%
  st_as_sf(crs = 4326, coords = c("lon", "lat")) %>%
  st_transform(crs_utm9)

bbox <- st_bbox(bc_coast)
coord <- coord_sf(xlim = c(bbox[1] - 6e4, bbox[3] - 3e5),
  ylim = c(bbox[2], bbox[4] - 5e5))

ggplot(bc_coast) +
  geom_sf() +
  geom_sf(data = surv, size = 0.5) +
  coord

surv_utm_coords <- st_coordinates(surv)
dat$X1000 <- surv_utm_coords[,1] / 1000
dat$Y1000 <- surv_utm_coords[,2] / 1000

spde <- sdmTMB::make_spde(dat, xy_cols = c("X1000", "Y1000"), n_knots = 700, type = "cutoff_search")
plot(spde)
mesh <- spde$mesh
tl <- length(mesh$graph$tv[, 1]) # the number of triangles in the mesh
posTri <- matrix(0, tl, 2)
for (i in seq_len(tl)) {
  temp <- mesh$loc[mesh$graph$tv[i, ], ]
  posTri[i, ] <- colMeans(temp)[c(1, 2)]
}
mesh_sf <- as.data.frame(posTri) %>%
  mutate(X = V1 * 1000, Y = V2 * 1000) %>%
  st_as_sf(crs = crs_utm9, coords = c("X", "Y"))

intersected <- st_intersects(mesh_sf, bc_coast)

mesh_df_water <- mesh_sf[lengths(intersected) == 0, ]
mesh_df_land <- mesh_sf[lengths(intersected) > 0, ]

pal <- RColorBrewer::brewer.pal(4, "Dark2")

ggplot(bc_coast) +
  geom_sf() +
  geom_sf(data = mesh_df_water, size = 1, colour = pal[3]) +
  geom_sf(data = mesh_df_land, size = 1, colour = pal[1]) +
  coord

water.triangles <- which(lengths(intersected) == 0)
land.triangles <- which(lengths(intersected) > 0)

plot(posTri[water.triangles,])
plot(posTri[land.triangles,])

barrier_spde <- INLA::inla.barrier.fem(mesh,
  barrier.triangles = land.triangles)

bspde <- spde
bspde$spde_barrier <- barrier_spde
dat$density1000 <- dat$density * 1000

mb <- sdmTMB(density1000 ~ 0 + as.factor(year),
  data = dat, spde = bspde, silent = FALSE, barrier = TRUE,
  barrier_scaling = c(1, 0.1))

mb

grid <- gfplot::synoptic_grid

original_time <- sort(unique(dat$year))
nd <- do.call("rbind",
  replicate(length(original_time), grid, simplify = FALSE))
nd[["year"]] <- rep(original_time, each = nrow(grid))

nd$X1000 <- nd$X
nd$Y1000 <- nd$Y

pb <- predict(mb, newdata = nd)

filter(pb, year %in% c(2017, 2018)) %>%
  ggplot(aes(X1000, Y1000, colour = est, fill = est)) +
  geom_tile(width = 2.1, height = 2.1) +
  facet_wrap(~year) +
  geom_point(data = filter(dat, year %in% c(2017, 2018)), aes(x = X1000, y = Y1000), inherit.aes = FALSE, pch = 4, size = 0.7, alpha = 0.4) +
  scale_colour_viridis_c() +
  scale_fill_viridis_c() +
  coord_fixed()

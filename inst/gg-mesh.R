library(ggplot2)

map_data <- rnaturalearth::ne_countries(
  scale = "large",
  returnclass = "sf", country = "canada"
)
# Crop the polygon for plotting and efficiency:
bc_coast <- suppressWarnings(suppressMessages(
  sf::st_crop(
    map_data,
    c(xmin = -134, ymin = 46, xmax = -120, ymax = 57)
  )
))
utm_zone9 <- 3156
bc_coast_proj <- sf::st_transform(bc_coast, crs = utm_zone9)
lims_x <- c(230957.7 + 105000, 1157991 - 570000) + c(-10000, 10000)
lims_y <- c(5366427 + 270000, 6353456 - 513000) + c(-10000, 10000)

# example mesh
dat <- sdmTMB::pcod # example data
dat$X1000 <- dat$X * 1000 # km to m
dat$Y1000 <- dat$Y * 1000 # km to m
mesh_1000 <- sdmTMB::make_mesh(dat, xy_cols = c("X1000", "Y1000"), cutoff = 5 * 1000)

g_mesh <- ggplot(bc_coast_proj) +
  geom_sf(colour = "grey86", lwd = 0.3, fill = "grey86") +
  inlabru::gg(mesh_1000$mesh, edge.color = "grey61") +
  coord_sf(xlim = lims_x, ylim = lims_y) +
  geom_point(data = dat, alpha = 0.8,
    mapping = aes(X1000, Y1000, size = density), colour = "#3182BD") +
  scale_size_area(max_size = 10) +
  theme_light() +
  theme(legend.position = "none") +
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank()
  )
g_mesh

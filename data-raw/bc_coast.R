# remotes::install_github("ropensci/rnaturalearthhires")
map_data <- rnaturalearth::ne_countries(
  scale = "large",
  returnclass = "sf", country = "canada")
bc_coast <- suppressWarnings(suppressMessages(
  st_crop(map_data,
    c(xmin = -134, ymin = 46, xmax = -120, ymax = 57))))

crs_utm9 <- 3156 # Pick a projection, here UTM9

sf::st_crs(bc_coast) <- 4326 # 'WGS84'; necessary on some installs

bc_coast <- sf::st_transform(bc_coast, crs_utm9)

# avoid Latin-1 NOTEs
# https://github.com/pbs-assess/sdmTMB/issues/158
bc_coast$name_es <- "Canada"
bc_coast$name_pt <- "Canada"

usethis::use_data(bc_coast, overwrite = TRUE)

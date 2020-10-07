# remotes::install_github("ropensci/rnaturalearthhires")
map_data <- rnaturalearth::ne_countries(
  scale = "large",
  returnclass = "sf", country = "canada")
bc_coast <- suppressWarnings(suppressMessages(
  st_crop(map_data,
    c(xmin = -134, ymin = 46, xmax = -120, ymax = 57))))
usethis::use_data(bc_coast)

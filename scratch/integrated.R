library(tinyVAST)
library(sf)
library(fmesher)
library(ggplot2)

# Load data
data(red_snapper)
data(red_snapper_shapefile)

# Plot data extent
plot(
  x = red_snapper$Lon,
  y = red_snapper$Lat,
  col = rainbow(3)[red_snapper$Data_type]
)
plot(red_snapper_shapefile, col = rgb(0, 0, 0, 0.2), add = TRUE)
legend("bottomleft",
  bty = "n",
  legend = levels(red_snapper$Data_type),
  fill = rainbow(3)
)

#
family <- list(
  "Encounter" = binomial(link = "cloglog"),
  "Count" = poisson(link = "log"),
  "Biomass_KG" = tweedie(link = "log")
)

# Relevel gear factor so Biomass_KG is base level
red_snapper$Data_type <- relevel(red_snapper$Data_type,
  ref = "Biomass_KG"
)

# Define mesh
mesh <- fm_mesh_2d(red_snapper[, c("Lon", "Lat")], cutoff = 0.5)

# define formula with a catchability covariate for gear
formula <- Response_variable ~ Data_type + factor(Year) + offset(log(AreaSwept_km2))

# make variable column
red_snapper$var <- "logdens"

# fit using tinyVAST
fit <- tinyVAST(
  data = red_snapper,
  formula = formula,
  sem = "logdens <-> logdens, sd_space",
  dsem = "logdens <-> logdens, 0, sd_spacetime",
  space_columns = c("Lon", "Lat"),
  spatial_graph = mesh,
  time_column = "Year",
  distribution_column = "Data_type",
  family = family,
  variable_column = "var",
  control = tinyVASTcontrol(silent = FALSE)
)

# now do it in sdmTMB

head(red_snapper)



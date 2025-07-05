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
.family <- list(
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
  space_term = "logdens <-> logdens, sd_space",
  spacetime_term = "logdens <-> logdens, 0, sd_spacetime",
  space_columns = c("Lon", "Lat"),
  spatial_domain = mesh,
  time_column = "Year",
  distribution_column = "Data_type",
  family = .family,
  variable_column = "var"
)

# now do it in sdmTMB

head(red_snapper)
distribution_column <- "Data_type"
data <- red_snapper

# match(data[, distribution_column], names(.family))

f <- .family[data[[distribution_column]]]
.valid_family[unlist(lapply(f, `[[`, "family"))]
.valid_link[unlist(lapply(f, `[[`, "link"))]

f <- list(poisson(), poisson(), delta_gamma(), delta_gamma(type = "poisson-link"))
ff <- lapply(f, `[[`, "family")
ff <- lapply(ff, \(x) if (length(x) == 1) c(x, NA) else x)
ff <- do.call(rbind, ff)

.valid_family[ff[, 1]]
.valid_family[ff[, 2]]
as.integer(!is.na(ff[, 2]))
as.integer(lapply(f, `[[`, "type") == "standard")
as.integer(lapply(f, `[[`, "type") == "poisson_link_delta")

.valid_link[family$link]

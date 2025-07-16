library(ggplot2)
library(fmesher)
library(sp)
library(showtext)

# Add Google Font
font_add_google("Roboto", "roboto")
font_add_google("Poppins", "poppins")
font_add_google("Roboto Condensed", "roboto condensed")
showtext_auto()

# Create hexagon vertices
center_x <- 0.5
center_y <- 0.5
radius <- 0.5

# Calculate the 6 vertices of the hexagon (rotated to have flat sides)
angles <- seq(0, 2*pi, length.out = 7) + pi/6  # 7 points to close the polygon, rotated 30 degrees
hex_x <- center_x + radius * cos(angles)
hex_y <- center_y + radius * sin(angles)

# Create hexagon data frame
hexagon_data <- data.frame(
  x = hex_x,
  y = hex_y
)

# Function to check if a point is inside the hexagon
point_in_hexagon <- function(x, y) {
  # Create a polygon object for the hexagon
  hex_poly <- Polygon(cbind(hex_x[-7], hex_y[-7]))  # Remove duplicate last point
  hex_poly_sp <- Polygons(list(hex_poly), ID = "hex")
  hex_spatial <- SpatialPolygons(list(hex_poly_sp))

  # Check if points are inside
  points_sp <- SpatialPoints(cbind(x, y))
  inside <- !is.na(over(points_sp, hex_spatial))
  return(inside)
}

# Generate random points inside the hexagon
# seed <- sample(1:1e6, 1)
set.seed(42)  # For reproducible results
n_points <- 100
attempts <- 0
max_attempts <- 1000

# Generate points inside hexagon
points_x <- numeric(0)
points_y <- numeric(0)

while(length(points_x) < n_points && attempts < max_attempts) {
  # Generate random points in the bounding box
  x_candidates <- runif(n_points * 2, 0, 1)
  y_candidates <- runif(n_points * 2, 0, 1)

  # Check which points are inside the hexagon
  inside_mask <- point_in_hexagon(x_candidates, y_candidates)

  # Add valid points
  valid_x <- x_candidates[inside_mask]
  valid_y <- y_candidates[inside_mask]

  points_x <- c(points_x, valid_x)
  points_y <- c(points_y, valid_y)

  attempts <- attempts + 1
}

# Take only the first n_points
points_x <- points_x[1:min(n_points, length(points_x))]
points_y <- points_y[1:min(n_points, length(points_y))]

obs <- data.frame(x = points_x, y = points_y, value = rlnorm(n_points))

# Add hexagon vertices to ensure triangulation covers the boundary
points_x <- c(points_x, hex_x[-7])  # Remove duplicate last point
points_y <- c(points_y, hex_y[-7])

# Create boundary for the hexagon using fmesher format
boundary <- fm_segm(cbind(hex_x, hex_y), is.bnd = TRUE)

# Create mesh using fmesher
mesh <- fm_mesh_2d(
  boundary = boundary,
  max.edge = 0.34,      # Maximum edge length
  cutoff = 0.03,       # Minimum distance between points
  offset = 0.11        # Offset for boundary
)

# Extract triangulation data from mesh
triangles <- mesh$graph$tv
points <- mesh$loc[, 1:2]

# Generate spatially correlated colors using Gaussian random field
seed <- sample(1:1e6, 1)
set.seed(463059)  # For reproducible colors
# set.seed(seed)  # For reproducible colors

# Create a Gaussian random field on the mesh vertices
n_vertices <- nrow(points)
# Simple approach: create spatial correlation based on distance
vertex_grf <- numeric(n_vertices)

# Add covariate effect: darker diagonally (60 degrees clockwise from horizontal)
# Rotate coordinates 60 degrees clockwise
angle <- -60 * pi / 180  # Convert to radians, negative for clockwise
rotation_matrix <- matrix(c(cos(angle), -sin(angle), sin(angle), cos(angle)), nrow = 2)

# Apply rotation to points
rotated_points <- t(rotation_matrix %*% t(points))

# Use the rotated x-coordinate for the covariate effect
x_rotated_normalized <- (rotated_points[, 1] - min(rotated_points[, 1])) / (max(rotated_points[, 1]) - min(rotated_points[, 1]))
covariate_effect <- 5 * x_rotated_normalized  # Adjust coefficient to control strength

# Generate correlated field using a simple exponential correlation
correlation_range <- 0.2  # Adjust this to control spatial correlation
for(i in 1:n_vertices) {
  # Base random value plus covariate effect
  vertex_grf[i] <- rnorm(1) + covariate_effect[i]

  # Add spatial correlation by averaging with nearby vertices
  distances <- sqrt(rowSums((points - matrix(points[i,], nrow = n_vertices, ncol = 2, byrow = TRUE))^2))
  weights <- exp(-distances / correlation_range)
  weights[i] <- 0  # Don't include self

  if(sum(weights) > 0) {
    # Add weighted contribution from nearby vertices (if any have been assigned)
    nearby_values <- vertex_grf[1:max(1, i-1)]
    nearby_weights <- weights[1:max(1, i-1)]
    if(length(nearby_values) > 0 && sum(nearby_weights) > 0) {
      vertex_grf[i] <- vertex_grf[i] + 0.5 * sum(nearby_values * nearby_weights) / sum(nearby_weights)
    }
  }
}

# Assign triangle colors as the mean of their vertex values
triangle_colors <- numeric(nrow(triangles))
for(i in 1:nrow(triangles)) {
  tri_indices <- triangles[i, ]
  triangle_colors[i] <- mean(vertex_grf[tri_indices])
}

# Create triangle data frame for ggplot with random colors
triangle_data <- data.frame()
for(i in 1:nrow(triangles)) {
  tri_indices <- triangles[i, ]
  tri_coords <- points[tri_indices, ]

  triangle_df <- data.frame(
    x = c(tri_coords[, 1], tri_coords[1, 1]),  # Close the triangle
    y = c(tri_coords[, 2], tri_coords[1, 2]),
    triangle_id = i,
    color_value = triangle_colors[i]
  )
  triangle_data <- rbind(triangle_data, triangle_df)
}

# Create the plot
border_col <- "grey25"
border_col <- RColorBrewer::brewer.pal(8, "Blues")[8]
viridis_option <- "G"
border_col <- viridis::viridis(option = viridis_option, n = 8)[2]

border_col <- "grey25"
border_col <- "#9e3434"

g <- ggplot() +
  geom_polygon(data = hexagon_data, aes(x = x, y = y),
    fill = NA, color = border_col, linewidth = 10) +
  # geom_polygon(data = triangle_data,
  #   aes(x = x, y = y, group = triangle_id, fill = color_value),
  #   color = "white", linewidth = 0) +
  # geom_polygon(data = triangle_data,
  #   aes(x = x, y = y, group = triangle_id), fill = "#FFFFFF20",
  #   color = "white", linewidth = 0) +
  # geom_polygon(data = triangle_data,
  #   aes(x = x, y = y, group = triangle_id), fill = NA,
  #   color = "#FFFFFF70", linewidth = 0.5) +
  geom_polygon(data = triangle_data,
    aes(x = x, y = y, group = triangle_id),
    color = "white", linewidth = 0.3, fill = "#D54545") +
  geom_polygon(data = hexagon_data, aes(x = x, y = y),
    fill = NA, color = border_col, linewidth = 0.3) +
  # geom_text(aes(x = 0.35, y = 0.8), label = "sdmTMB",
            # size = 13, color = "white", angle = 30, family = "roboto condensed") +

  # geom_text(aes(x = 0.5, y = 0.58), label = "sdmTMB",
  # size = 14, color = "grey25", angle = 0, family = "roboto condensed") +

  geom_text(aes(x = 0.47, y = 0.57), label = "sdmTMB",
    size = 16, color = "white", angle = 30, family = "roboto condensed") +

  # geom_text(aes(x = 0.343, y = 0.755), label = "sdmTMB",
    # size = 16, color = "white", angle = 30, family = "roboto condensed", fontface = "italic") +

    # geom_text(aes(x = 0.5, y = 0.56), label = "sdmTMB",
    # size = 13, color = "grey25", angle = 0, family = "roboto condensed") +

  # Filter out points that overlap with the text area
  # geom_point(data = obs[!(obs$x < 0.65 & obs$y > 0.55), ], aes(x, y, size = value), colour = "white", alpha = 0.8, pch = 20) +
  # Filter out points that overlap with the diagonal text area
  # Create a diagonal exclusion zone around the text "sdmTMB"
  # Create a wider diagonal exclusion zone around the text "sdmTMB" at 30 degrees
  geom_point(data = obs[!(abs((obs$y - 0.57) - tan(30 * pi/180) * (obs$x - 0.47)) < 0.12 & 
                          obs$x > 0.2 & obs$x < 0.8 & obs$y > 0.35 & obs$y < 0.8), ], 
             aes(x, y, size = value), colour = "white", alpha = 0.8, pch = 20) +
  # scale_fill_viridis_c(option = viridis_option, direction = -1) +
  # scale_fill_distiller(palette = "Reds", direction = 1, ) +
  # scale_fill_gradient(low = "white", high = RColorBrewer::brewer.pal(8, "Blues")[8]) +
  scale_size_area(max_size = 4) +
  guides(size = FALSE, fill = FALSE) +
  coord_fixed(ratio = 1) +
  theme_minimal()
  # theme_void()
print(g)


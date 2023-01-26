#' Transform a mesh object into a mesh with correlation barriers
#'
#' @param spde_obj Output from [make_mesh()].
#' @param barrier_sf An sf object with polygons defining the barriers. For
#'   example, a coastline dataset for ocean data. **Note that this object must
#'   have the same projection as the data used to generate the x and y columns
#'   in `spde_obj`.**
#' @param range_fraction The fraction of the spatial range that barrier
#'   triangles have.
#' @param proj_scaling If `spde_obj` was created with scaling of the coordinates
#'   after the projection (e.g., dividing UTMs by 1000 so the spatial range is
#'   on a reasonable scale) the x and y values in `spde_obj` are multiplied by
#'   this scaling factor before applying the projection from `barrier_sf`.
#' @param plot Logical.
#'
#' @return A list similar to [make_mesh()] but with `spde_barrier` and a
#' couple other helper list elements added.
#'
#' If `plot = TRUE`, then a basic plot will be created as a side effect. Each
#' grey dot represents the center of a "normal" mesh triangle. Each red cross
#' represents the center of a "barrier" mesh triangle.
#' @export
#' @references
#' Bakka, H., Vanhatalo, J., Illian, J., Simpson, D., and Rue, H. 2019.
#' Non-stationary Gaussian models with physical barriers.
#' <https://arxiv.org/abs/1608.03787>
#'
#' <https://sites.google.com/a/r-inla.org/www/barrier-model>
#'
#' <https://haakonbakkagit.github.io/btopic107.html>
#' @examples
#' if (require("sf", quietly = TRUE) &&
#'   require("ggplot2", quietly = TRUE) &&
#'   require("dplyr", quietly = TRUE) &&
#'   require("INLA", quietly = TRUE)) {
#'
#' # First, download coastline data for our region.
#' # We will use `bc_coast` from the package data,
#' # but you can recreate it with the following.
#'
#' # For applied situations on finer scales, you may with to use scale = "large".
#' # For that, first: remotes::install_github("ropensci/rnaturalearthhires")
#' # map_data <- rnaturalearth::ne_countries(
#' #   scale = "medium",
#' #   returnclass = "sf", country = "canada")
#' #
#' # # Crop the polygon for plotting and efficiency:
#' # st_bbox(map_data)
#' # bc_coast <- suppressWarnings(suppressMessages(
#' #   st_crop(map_data,
#' #     c(xmin = -134, ymin = 46, xmax = -120, ymax = 57))))
#'
#' crs_utm9 <- 3156 # Pick a projection, here UTM9
#'
#' st_crs(bc_coast) <- 4326 # 'WGS84'; necessary on some installs
#' bc_coast <- st_transform(bc_coast, crs_utm9)
#'
#' # Project our survey data coordinates:
#' survey <- pcod %>% select(lon, lat, density) %>%
#'   st_as_sf(crs = 4326, coords = c("lon", "lat")) %>%
#'   st_transform(crs_utm9)
#'
#' # Plot our coast and survey data:
#' ggplot(bc_coast) +
#'   geom_sf() +
#'   geom_sf(data = survey, size = 0.5)
#'
#' # Note that a barrier mesh won't do much here for this
#' # example data set, but we nonetheless use it as an example.
#'
#' # Prepare for making the mesh
#' # First, we will extract the coordinates:
#' surv_utm_coords <- st_coordinates(survey)
#'
#' # Then we will scale coordinates to km so the range parameter
#' # is on a reasonable scale for estimation:
#' pcod$X1000 <- surv_utm_coords[,1] / 1000
#' pcod$Y1000 <- surv_utm_coords[,2] / 1000
#'
# # Construct our mesh:
#' spde <- make_mesh(pcod, xy_cols = c("X1000", "Y1000"),
#'   n_knots = 200, type = "kmeans")
#' plot(spde)
#'
#' # Add on the barrier mesh component:
#' bspde <- add_barrier_mesh(
#'   spde, bc_coast, range_fraction = 0.1,
#'   proj_scaling = 1000, plot = TRUE
#' )
#'
#' # In the above, the grey dots are the centre of triangles that are in the
#' # ocean. The red crosses are centres of triangles that are over land. The
#' # spatial range will be assumed to be 0.1 (`range_fraction`) over land compared
#' # to over water.
#'
#' # We can make a more advanced plot if we want:
#' mesh_df_water <- bspde$mesh_sf[bspde$normal_triangles, ]
#' mesh_df_land <- bspde$mesh_sf[bspde$barrier_triangles, ]
#' ggplot(bc_coast) +
#'   geom_sf() +
#'   geom_sf(data = mesh_df_water, size = 1, colour = "blue") +
#'   geom_sf(data = mesh_df_land, size = 1, colour = "green")
#'
#' # Now, when we fit our model with the new mesh, it will automatically
#' # include a barrier structure in the spatial correlation:
#' fit <- sdmTMB(density ~ s(depth, k = 3), data = pcod, mesh = bspde,
#'   family = tweedie(link = "log"))
#' fit
#' }

add_barrier_mesh <- function(spde_obj, barrier_sf, range_fraction = 0.2,
  proj_scaling = 1, plot = FALSE) {

  if (!requireNamespace("INLA", quietly = TRUE)) {
    cli_abort("INLA must be installed to use this function.")
  }
  if (!requireNamespace("sf", quietly = TRUE)) {
    cli_abort("The sf package must be installed to use this function.")
  }

  assert_that(
    is.numeric(range_fraction), range_fraction <= 1,
    range_fraction > 0, length(range_fraction) == 1
  )
  assert_that(is.numeric(proj_scaling), length(proj_scaling) == 1)
  assert_that("sf" %in% class(barrier_sf))
  assert_that("sdmTMBmesh" %in% class(spde_obj))
  assert_that(is.logical(plot))

  mesh <- spde_obj$mesh
  tl <- length(mesh$graph$tv[, 1]) # the number of triangles in the mesh
  pos_tri <- matrix(0, tl, 2)
  for (i in seq_len(tl)) {
    temp <- mesh$loc[mesh$graph$tv[i, ], ]
    pos_tri[i, ] <- colMeans(temp)[c(1, 2)]
  }
  mesh_sf <- as.data.frame(pos_tri)
  mesh_sf$X <- mesh_sf$V1 * proj_scaling
  mesh_sf$Y <- mesh_sf$V2 * proj_scaling
  mesh_sf <- sf::st_as_sf(mesh_sf, crs = sf::st_crs(barrier_sf), coords = c("X", "Y"))

  intersected <- sf::st_intersects(mesh_sf, barrier_sf)
  water.triangles <- which(lengths(intersected) == 0)
  land.triangles <- which(lengths(intersected) > 0)
  # mesh_df_water <- mesh_sf[water.triangles, ]
  # mesh_df_land <- mesh_sf[land.triangles, ]

  if (plot) {
    plot(pos_tri[water.triangles, ],
      col = "grey40", asp = 1,
      xlab = spde_obj$xy_cols[1], ylab = spde_obj$xy_cols[2]
    )
    points(pos_tri[land.triangles, ], col = "red", pch = 4)
    # g <- ggplot2::ggplot(barrier_sf) +
    #   ggplot2::geom_sf() +
    #   ggplot2::geom_sf(data = mesh_df_water, size = 1, colour = "blue") +
    #   ggplot2::geom_sf(data = mesh_df_land, size = 1, colour = "green")
    # print(g)
  }
  barrier_spde <- INLA::inla.barrier.fem(mesh, barrier.triangles = land.triangles)
  spde_obj$spde_barrier <- barrier_spde
  spde_obj$barrier_scaling <- c(1, range_fraction)
  spde_obj$mesh_sf <- mesh_sf
  spde_obj$barrier_triangles <- land.triangles
  spde_obj$normal_triangles <- water.triangles
  spde_obj
}

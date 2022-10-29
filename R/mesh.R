#' Construct an SPDE mesh for sdmTMB
#'
#' Construct an SPDE mesh for use with sdmTMB.
#'
#' @param data A data frame.
#' @param xy_cols A character vector of x and y column names contained in
#'   `data`. These should likely be in an equal distance projection. For
#'   a helper function to convert to UTMs, see [add_utm_columns()].
#' @param type Method to create the mesh. Also see `mesh` argument to supply
#'   your own mesh.
#' @param cutoff An optional cutoff if type is `"cutoff"`. "The minimum allowed
#'   distance between points in the mesh". See [INLA::inla.mesh.create()].
#'   Smaller values create meshes with more knots. Points further apart than this
#'   value will receive a separate vertex in the mesh before any mesh refinement.
#' @param n_knots The number of desired knots if `type` is not `"cutoff"`.
#' @param seed Random seed. Affects [stats::kmeans()] determination of knot
#'   locations if `type = "kmeans"`.
#' @param refine Logical or list to pass to [INLA::inla.mesh.create()].
#' @param mesh An optional mesh created via INLA instead of using the above
#'   convenience options.
#'
#' @return
#' `make_mesh()`: A list of class `sdmTMBmesh`. The element `mesh` is the output from
#' [INLA::inla.mesh.create()] and the element `spde` is the output from
#' [INLA::inla.spde2.matern()].
#'
#' @export
#'
#' @examplesIf inla_installed()
#' mesh <- make_mesh(pcod, c("X", "Y"), cutoff = 30, type = "cutoff")
#' plot(mesh)
#'
#' \donttest{
#' mesh <- make_mesh(pcod, c("X", "Y"), cutoff = 5, type = "cutoff")
#' plot(mesh)
#'
#' mesh <- make_mesh(pcod, c("X", "Y"), n_knots = 50, type = "cutoff_search")
#' plot(mesh)
#'
#' mesh <- make_mesh(pcod, c("X", "Y"), n_knots = 50, type = "kmeans")
#' plot(mesh)
#'
#' # Defining a mesh directly with INLA:
#' bnd <- INLA::inla.nonconvex.hull(cbind(pcod$X, pcod$Y), convex = -0.05)
#' inla_mesh <- INLA::inla.mesh.2d(
#'   boundary = bnd,
#'   max.edge = c(20, 50),
#'   offset = -0.05,
#'   cutoff = c(2, 5),
#'   min.angle = 10
#' )
#' mesh <- make_mesh(pcod, c("X", "Y"), mesh = inla_mesh)
#' plot(mesh)
#' }
make_mesh <- function(data, xy_cols,
                      type = c("kmeans", "cutoff", "cutoff_search"),
                      cutoff, n_knots,
                      seed = 42,
                      refine = list(min.angle = 21, max.edge = Inf, max.n.strict = -1, max.n = 1000),
                      mesh = NULL) {

  if (!requireNamespace("INLA", quietly = TRUE)) {
    cli_abort("INLA must be installed to use this function.")
  }

  if (missing(xy_cols) || is.numeric(xy_cols) || is.numeric(data)) {
    msg <- paste0(
      "It looks like you are using an old format of make_mesh(). ",
      "The function now uses `data` and `xy_cols` arguments ",
      "to enable carrying through the x and y column names ",
      "to the predict function. Please update your code."
    )
    cli_abort(msg)
  }

  if (max(data[[xy_cols[1]]]) > 1e4 || max(data[[xy_cols[2]]] > 1e4)) {
    msg <- paste0(
      "The x or y column values are fairly large. ",
      "This can cause estimation problems since the spatial range ",
      "is dependent on the scale of the coordinates. ",
      "Consider scaling the x and y coordinates. ",
      "For example, try working in UTM km instead of UTM m by divided by 1000."
      )
    cli_warn(msg)
  }
  type <- match.arg(type)
  if (!missing(n_knots) && type == "cutoff") {
    msg <- paste0(
      "`n_knots` specified but `type = 'cutoff'. ",
      "The default mesh type now uses a `cutoff` ",
      "intead of the number of knots. You can obtain the previous default ",
      "mesh type with `type = 'kmeans'`."
      )
    cli_abort(msg)
  }

  if (!missing(cutoff) && missing(n_knots)) {
    type <- "cutoff"
  }
  if (missing(cutoff) && type == "cutoff" && is.null(mesh)) {
    cli_abort("You need to specify the `cutoff` argument.")
  }
  if (missing(n_knots) && type != "cutoff" && is.null(mesh)) {
    cli_abort("You need to specify the `n_knots` argument.")
  }
  loc_xy <- as.matrix(data[, xy_cols, drop = FALSE])
  loc_centers <- NA

  if (is.null(mesh)) {
    if (type == "kmeans") {
      if (n_knots >= nrow(loc_xy)) {
        cli_warn(
          "Reducing `n_knots` to be one less than the number of data points."
        )
        n_knots <- nrow(loc_xy) - 1
      }
      set.seed(seed)
      knots <- stats::kmeans(x = loc_xy, centers = n_knots)
      loc_centers <- knots$centers
      mesh <- INLA::inla.mesh.create(loc_centers, refine = TRUE)
    } else if (type == "cutoff") {
      mesh <- INLA::inla.mesh.create(loc_xy, refine = TRUE, cutoff = cutoff)
    } else {
      mesh <- binary_search_knots(loc_xy, n_knots = n_knots, refine = refine)
    }
  } else {
    knots <- list()
    loc_centers <- NA
  }
  spde <- INLA::inla.spde2.matern(mesh)
  A <- INLA::inla.spde.make.A(mesh, loc = loc_xy)

  fake_data <- data
  fake_data[["sdm_spatial_id"]] <- seq(1, nrow(fake_data))

  structure(list(
    loc_xy = loc_xy, xy_cols = xy_cols, mesh = mesh, spde = spde,
    loc_centers = loc_centers, A_st = A,
    sdm_spatial_id = fake_data$sdm_spatial_id
  ), class = "sdmTMBmesh")}

binary_search_knots <- function(loc_xy,
                                n_knots,
                                min = 1e-4,
                                max = 1e4,
                                length = 1e6,
                                refine = TRUE) {
  vec <- exp(seq(log(min), log(max), length.out = length))
  L <- 0
  R <- length(vec)
  realized_knots <- Inf
  while (L <= R) {
    m <- floor((L + R) / 2)
    mesh <- INLA::inla.mesh.create(loc_xy, refine = refine, cutoff = vec[m])
    realized_knots <- mesh$n
    pretty_cutoff <- sprintf("%.2f", round(vec[m], 2))
    cat("cutoff =", pretty_cutoff, "| knots =", realized_knots)
    if (realized_knots > n_knots) {
      L <- m + 1
      cat(" |", clisymbols::symbol$arrow_down, "\n")
    } else if (realized_knots < n_knots) {
      R <- m - 1
      cat(" |", clisymbols::symbol$arrow_up, "\n")
    } else {
      cat(" |", clisymbols::symbol$tick, "\n")
      return(mesh)
    }
  }
  cat("cutoff =", pretty_cutoff, "| knots =", mesh$n, "\n")
  mesh
}

#' @param x Output from [make_mesh()].
#' @param ... Passed to [graphics::plot()].
#'
#' @importFrom graphics points
#' @return `plot.sdmTMB()`: A plot of the mesh and data points.
#' @rdname make_mesh
#' @export
plot.sdmTMBmesh <- function(x, ...) {
  plot(x$mesh, main = NA, edge.color = "grey60", asp = 1, ...)
  points(x$loc_xy, pch = 21, col = "#00000070")
  points(x$loc_centers, pch = 20, col = "red")
}

# from TMB examples repository:
make_anisotropy_spde <- function(spde, anistropy = TRUE) {
  if (anistropy) {
    inla_mesh <- spde$mesh
    Dset <- 1:2
    inla_mesh <- spde$mesh
    # Triangle info
    TV <- inla_mesh$graph$tv # Triangle to vertex indexing
    V0 <- inla_mesh$loc[TV[, 1], Dset] # V = vertices for each triangle
    V1 <- inla_mesh$loc[TV[, 2], Dset]
    V2 <- inla_mesh$loc[TV[, 3], Dset]
    E0 <- V2 - V1 # E = edge for each triangle
    E1 <- V0 - V2
    E2 <- V1 - V0
    # Calculate Areas
    TmpFn <- function(Vec1, Vec2) abs(det(rbind(Vec1, Vec2)))
    Tri_Area <- rep(NA, nrow(E0))
    for (i in seq_len(length(Tri_Area))) Tri_Area[i] <- TmpFn(E0[i, ], E1[i, ]) / 2

    ret <- list(
      n_s = spde$spde$n.spde,
      n_tri = nrow(TV),
      Tri_Area = Tri_Area,
      E0 = E0,
      E1 = E1,
      E2 = E2,
      TV = TV - 1,
      G0 = spde$spde$param.inla$M0,
      G0_inv = as(Matrix::diag(1 / Matrix::diag(spde$spde$param.inla$M0)), "TsparseMatrix")
    )
  } else {
    ret <- list(
      n_s = 0L, n_tri = 0L, Tri_Area = rep(0, 1), E0 = matrix(0, 1),
      E1 = matrix(0, 1), E2 = matrix(0, 1), TV = matrix(0, 1),
      G0 = Matrix::Matrix(0, 1, 1, doDiag = FALSE), G0_inv = Matrix::Matrix(0, 1, 1, doDiag = FALSE)
    )
  }
  ret
}

make_barrier_spde <- function(spde) {
  if ("spde_barrier" %in% names(spde)) {
    C0 <- spde$spde_barrier$C[[1]]
    C1 <- spde$spde_barrier$C[[2]]
    D0 <- spde$spde_barrier$D[[1]]
    D1 <- spde$spde_barrier$D[[2]]
    .I <- spde$spde_barrier$I
  } else {
    C0 <- rep(1, 2)
    C1 <- rep(1, 2)
    D0 <- Matrix::Matrix(0, 1, 1, doDiag = FALSE)
    D1 <- Matrix::Matrix(0, 1, 1, doDiag = FALSE)
    .I <- Matrix::Matrix(0, 1, 1, doDiag = FALSE)
  }
  list(C0 = C0, C1 = C1, D0 = D0, D1 = D1, I = .I)
}


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

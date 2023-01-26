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


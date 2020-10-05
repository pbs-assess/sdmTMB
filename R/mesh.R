#' Construct an SPDE mesh
#'
#' @param x X numeric vector.
#' @param y Y numeric vector.
#' @param n_knots The number of knots.
#' @param seed Random seed. Affects [stats::kmeans()] determination of knot locations.
#' @param mesh An optional mesh created via INLA. If supplied, this mesh will be
#'   used instead of creating one with [stats::kmeans()] and the `n_knots`
#'   argument.
#'
#' @importFrom graphics points
#' @export
#' @examples
#' sp <- make_spde(pcod$X, pcod$Y, n_knots = 25)
#' plot_spde(sp)
#' \donttest{
#' loc_xy <- cbind(pcod$X, pcod$Y)
#' bnd <- INLA::inla.nonconvex.hull(as.matrix(loc_xy), convex = -0.05)
#' mesh <- INLA::inla.mesh.2d(
#'   boundary = bnd,
#'   max.edge = c(20, 50),
#'   offset = -0.05,
#'   cutoff = c(2, 5),
#'   min.angle = 10
#' )
#' sp2 <- make_spde(pcod$X, pcod$Y, mesh = mesh)
#' plot_spde(sp2)
#' }
#'
#' # make_spde <- function(x, y, n_knots, seed = 42, mesh = NULL) {
#' #   loc_xy <- cbind(x, y)
#' #
#' #   if (is.null(mesh)) {
#' #     if (n_knots >= nrow(loc_xy)) {
#' #       warning(
#' #         "Reducing `n_knots` to be one less than the ",
#' #         "number of data points."
#' #       )
#' #       n_knots <- nrow(loc_xy) - 1
#' #     }
#' #     set.seed(seed)
#' #     knots <- stats::kmeans(x = loc_xy, centers = n_knots)
#' #     loc_centers <- knots$centers
#' #     mesh <- INLA::inla.mesh.create(loc_centers, refine = TRUE)
#' #   } else {
#' #     knots <- list()
#' #     loc_centers <- NA
#' #   }
#' #   spde <- INLA::inla.spde2.matern(mesh)
#' #   A <- INLA::inla.spde.make.A(mesh, loc = loc_xy)
#' #   list(
#' #     x = x, y = y, mesh = mesh, spde = spde,
#' #     loc_centers = loc_centers, A = A
#' #   )
#' # }


#' Construct an SPDE mesh
#'
#' @param data A data frame.
#' @param xy_cols A character vector of x and y column names contained in
#'   `data`.
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
#' @return A list of class sdmTMBmesh.
#' @export
#'
#' @examples
#' sp <- make_spde(pcod, c("X", "Y"), n_knots = 50, type = "cutoff-search")
#' plot(sp)
#'
#' sp <- make_spde(pcod, c("X", "Y"), n_knots = 50, type = "kmeans")
#' plot(sp)
#'
#' sp <- make_spde(pcod, c("X", "Y"), cutoff = 30, type = "cutoff")
#' plot(sp)
#'
#' sp <- make_spde(pcod, c("X", "Y"), cutoff = 5, type = "cutoff")
#' plot(sp)
#'
#' \donttest{
#' # Defining a mesh directly with INLA:
#' loc_xy <- cbind(pcod$X, pcod$Y)
#' bnd <- INLA::inla.nonconvex.hull(loc_xy, convex = -0.05)
#' mesh <- INLA::inla.mesh.2d(
#'   boundary = bnd,
#'   max.edge = c(20, 50),
#'   offset = -0.05,
#'   cutoff = c(2, 5),
#'   min.angle = 10
#' )
#' sp2 <- make_spde(pcod, c("X", "Y"), mesh = mesh)
#' plot(sp2)
#' }
#'
make_spde <- function(data, xy_cols,
                      type = c("cutoff", "cutoff-search", "kmeans"),
                      cutoff, n_knots,
                      seed = 42,
                      refine = list(min.angle = 21, max.edge = Inf, max.n.strict = -1, max.n = 1000),
                      mesh = NULL) {
  if (missing(xy_cols) || is.numeric(xy_cols) || is.numeric(data)) {
    stop("It looks like you are using an old format of make_spde(). ",
      "The function now uses `data` and `xy_cols` arguments ",
      "to enable carrying through the x and y column names ",
      "to the predict function. Please update your code. ",
      "Also note that the default mesh type now uses a `cutoff` ",
      "intead of the number of knots. You can obtain the previous default ",
      "mesh type with `type = 'kmeans'`.",
      call. = FALSE
    )
  }
  type <- match.arg(type)
  if (!missing(n_knots) && type == "cutoff") {
    stop("`n_knots` specified but `type = 'cutoff'. ",
      "The default mesh type now uses a `cutoff` ",
      "intead of the number of knots. You can obtain the previous default ",
      "mesh type with `type = 'kmeans'`.", call. = FALSE)
  }
  if (missing(cutoff) && type == "cutoff") {
    stop("You need to specify the `cutoff` argument.", call. = FALSE)
  }
  if (missing(n_knots) && type != "cutoff") {
    stop("You need to specify the `n_knots` argument.", call. = FALSE)
  }
  loc_xy <- as.matrix(data[, xy_cols, drop = FALSE])
  loc_centers <- NA

  if (is.null(mesh)) {
    if (type == "kmeans") {
      if (n_knots >= nrow(loc_xy)) {
        warning(
          "Reducing `n_knots` to be one less than the number of data points.",
          call. = FALSE
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
  structure(list(
    loc_xy = loc_xy, xy_cols = xy_cols, mesh = mesh, spde = spde,
    loc_centers = loc_centers, A = A
  ), class = "sdmTMBmesh")
}

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
  cat("cutoff =", pretty_cutoff, "| knots =", mesh$n, "¯\\_(ツ)_/¯\n")
  mesh
}


#' @param object Output from [make_spde()].
#' @rdname make_spde
#' @importFrom graphics plot
#' @export
plot_spde <- function(object) {
  .Deprecated("plot")
  plot(object$mesh, main = NA, edge.color = "grey60", asp = 1)
  if ("x" %in% names(object)) {
    points(object$x, object$y, pch = 21, col = "#00000070")
  } else {
    points(object$loc_xy, pch = 21, col = "#00000070")
  }
  points(object$loc_centers, pch = 20, col = "red")
}

#' Plot SPDE mesh object
#'
#' @param object Output from [make_spde2()].
#'
#' @return A plot
#' @export
plot.sdmTMBmesh <- function(object) {
  plot(object$mesh, main = NA, edge.color = "grey60", asp = 1)
  points(object$loc_xy, pch = 21, col = "#00000070")
  points(object$loc_centers, pch = 20, col = "red")
}

# from TMB examples repository:
make_anisotropy_spde <- function(spde) {
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
  for (i in 1:length(Tri_Area)) Tri_Area[i] <- TmpFn(E0[i, ], E1[i, ]) / 2

  list(
    n_s = spde$spde$n.spde,
    n_tri = nrow(TV),
    Tri_Area = Tri_Area,
    E0 = E0,
    E1 = E1,
    E2 = E2,
    TV = TV - 1,
    G0 = spde$spde$param.inla$M0,
    G0_inv = as(Matrix::diag(1 / Matrix::diag(spde$spde$param.inla$M0)), "dgTMatrix")
  )
}

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
#'
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

make_spde <- function(x, y, n_knots, seed = 42, mesh = NULL) {
  loc_xy <- cbind(x, y)

  if (is.null(mesh)) {
    if (n_knots >= nrow(loc_xy)) {
      warning(
        "Reducing `n_knots` to be one less than the ",
        "number of data points."
      )
      n_knots <- nrow(loc_xy) - 1
    }
    set.seed(seed)
    knots <- stats::kmeans(x = loc_xy, centers = n_knots)
    loc_centers <- knots$centers
    mesh <- INLA::inla.mesh.create(loc_centers, refine = TRUE)
  } else {
    knots <- list()
    knots$cluster <- vapply(seq_len(nrow(loc_xy)), function(i)
      RANN::nn2(mesh$loc[, 1:2, drop = FALSE],
        t(as.numeric(loc_xy[i, , drop = FALSE])),
        k = 1L
      )$nn.idx,
      FUN.VALUE = 1L
    )
    loc_centers <- NA
  }
  spde <- INLA::inla.spde2.matern(mesh)
  list(
    x = x, y = y, mesh = mesh, spde = spde, cluster = knots$cluster,
    loc_centers = loc_centers
  )
}

#' @param object Output from [make_spde()].
#' @rdname make_spde
#' @importFrom graphics plot
#' @export
plot_spde <- function(object) {
  plot(object$mesh, main = NA, edge.color = "grey60", asp = 1)
  points(object$x, object$y, pch = 21, col = "#00000070")
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
    n_s      = spde$spde$n.spde,
    n_tri    = nrow(TV),
    Tri_Area = Tri_Area,
    E0       = E0,
    E1       = E1,
    E2       = E2,
    TV       = TV - 1,
    G0       = spde$spde$param.inla$M0,
    G0_inv   = as(Matrix::diag(1/Matrix::diag(spde$spde$param.inla$M0)), "dgTMatrix"))
}

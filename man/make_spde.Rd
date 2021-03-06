% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mesh.R
\name{make_spde}
\alias{make_spde}
\alias{plot_spde}
\title{Construct an SPDE mesh (depreciated)}
\usage{
make_spde(x, y, n_knots, seed = 42, mesh = NULL)

plot_spde(object)
}
\arguments{
\item{x}{X numeric vector.}

\item{y}{Y numeric vector.}

\item{n_knots}{The number of knots.}

\item{seed}{Random seed. Affects \code{\link[stats:kmeans]{stats::kmeans()}} determination of knot locations.}

\item{mesh}{An optional mesh created via INLA. If supplied, this mesh will be
used instead of creating one with \code{\link[stats:kmeans]{stats::kmeans()}} and the \code{n_knots}
argument.
up}

\item{object}{Output from \code{\link[=make_spde]{make_spde()}}.}
}
\description{
\strong{Depreciated; please use \code{\link[=make_mesh]{make_mesh()}} instead.}
}
\examples{
# **Depreciated; please use `make_mesh()` instead.**
\donttest{
sp <- make_spde(pcod$X, pcod$Y, n_knots = 25)
plot_spde(sp)

loc_xy <- cbind(pcod$X, pcod$Y)
bnd <- INLA::inla.nonconvex.hull(as.matrix(loc_xy), convex = -0.05)
mesh <- INLA::inla.mesh.2d(
  boundary = bnd,
  max.edge = c(20, 50),
  offset = -0.05,
  cutoff = c(2, 5),
  min.angle = 10
)
sp2 <- make_spde(pcod$X, pcod$Y, mesh = mesh)
plot_spde(sp2)
}
}

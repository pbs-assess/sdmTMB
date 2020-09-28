# https://en.wikipedia.org/wiki/Binary_search_algorithm

binary_search <- function(x, vec) {
  L <- 0
  R <- length(vec)
  while (L <= R) {
    m <- floor((L + R) / 2)
    cat("m =", m, "\n")
    if (vec[m] < x) {
      L <- m + 1
    } else if (vec[m] > x) {
      R <- m - 1
    } else {
      return(m)
    }
  }
}

binary_search(500, seq(1, 10000))

#' Binary search for an INLA mesh with a cutoff that matches the desired knots
#'
#' @param x A vector of x coordinates.
#' @param y A vector of y coordinates.
#' @param n_knots The desired number of knots.
#' @param min The minimum cutoff evaluated.
#' @param max The maximum cutoff evaluated.
#' @param length The number of cutoffs evaluated between `min` and `max` where the increments are log distributed.
#'
#' @return A mesh from [INLA::inla.mesh.create()].

binary_search_knots <- function(x, y,
                                n_knots,
                                min = 1e-4,
                                max = 1e4,
                                length = 1e6) {
  vec <- exp(seq(log(min), log(max), length.out = length))
  loc_xy <- cbind(x, y)
  L <- 0
  R <- length(vec)
  realized_knots <- Inf
  while (L <= R) {
    m <- floor((L + R) / 2)
    mesh <- INLA::inla.mesh.create(loc_xy, refine = TRUE, cutoff = vec[m])
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

library(sdmTMB)
x <- pcod$X
y <- pcod$Y
s <- binary_search_knots(x, y, 500)
plot(s)
points(x, y, col = "#FF000090", pch = 20, cex = 0.5)
s <- binary_search_knots(x, y, 178)
s <- binary_search_knots(x, y, 200)
# an example that gets close but will fail:
s <- binary_search_knots(x, y, 500, length = 1e3)

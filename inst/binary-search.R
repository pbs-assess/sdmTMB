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

binary_search_knots <- function(loc_xy, n_knots, min = 1e-4, max = 1e4, length = 1e6) {
  vec <- exp(seq(log(min), log(max), length.out = length))
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

s <- binary_search_knots(loc_xy, 500)
plot(s)
points(loc_xy, col = "#FF000030")
s <- binary_search_knots(loc_xy, 178)
s <- binary_search_knots(loc_xy, 200)
s <- binary_search_knots(loc_xy, 500, length = 1e3)

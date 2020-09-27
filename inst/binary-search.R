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

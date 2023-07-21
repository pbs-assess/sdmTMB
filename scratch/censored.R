library(dplyr)
library(ggplot2)
theme_set(theme_light())

dcens_pois <- function(x, lambda) {
  1 - ppois(x - 1, lambda)
}

dcens_pois_upper <- function(x, lambda, upper) {
  ppois(upper, lambda) - ppois(x - 1, lambda)
}

# x <- seq(1, 40)
# mu <- 30
# u <- 35
#
# dc <- dcens_pois(x, 30)
# dcu <- dcens_pois_upper(x, 30, upper = u)
# d <- dpois(x, 30)
#
# dat <- data.frame(
#   x = c(x, x, x),
#   dens = c(d, dc, dcu),
#   type = c(
#     rep("dpois", length(x)),
#     rep("dcens_pois", length(x)),
#     rep("dcens_pois_upper", length(x))
#   )
# )
#
# ggplot(dat, aes(x, dens))+ facet_wrap(~type) + geom_col() +
#   geom_vline(xintercept = 15, col = "red", lwd = 2)

# ------------------------------------------------------------------

mu <- seq(1, 100)

df <- expand.grid(x = seq(0, 30, 5), u = c(seq(0, 60, 10), Inf)) |>
  dplyr::filter(u >= x)

out <- purrr::pmap_dfr(df, function(x, u) {
  dc <- dcens_pois(x, mu)
  dcu <- dcens_pois_upper(x, mu, upper = u)
  d <- dpois(x, mu)
  dat <- data.frame(
    x = paste("obs =", x),
    u = paste("upper =", u),
    xval = x,
    uval = u,
    mu = c(mu, mu, mu),
    dens = c(d, dc, dcu),
    type = c(
      rep("dpois", length(mu)),
      rep("dcens_pois", length(mu)),
      rep("dcens_pois_upper", length(mu))
    )
  )
  dat
})

dplyr::filter(out, type == "dpois") |>
  ggplot(aes(mu, dens))+ facet_grid(x~u) + geom_col(col = "grey50", fill = "grey50") +
  geom_vline(aes(xintercept = xval), col = "black", lty = 2)

dplyr::filter(out, type == "dcens_pois_upper") |>
  ggplot(aes(mu, dens))+ facet_grid(x~u) + geom_col(col = "grey50", fill = "grey50") +
  geom_vline(aes(xintercept = uval), col = "red") +
  geom_vline(aes(xintercept = xval), col = "black", lty = 2)

# dplyr::filter(out, type == "dcens_pois") |>
#   ggplot(aes(mu, dens))+ facet_grid(vars(x)) + geom_col(col = "grey50", fill = "grey50") +
#   geom_vline(aes(xintercept = xval), col = "black", lty = 2)

dplyr::filter(out, type == "dcens_pois_upper") |>
  ggplot(aes(mu, log(dens)))+ facet_grid(x~u) + geom_line() +
  geom_vline(aes(xintercept = uval), col = "red") +
  geom_vline(aes(xintercept = xval), col = "black", lty = 2)


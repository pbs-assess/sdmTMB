dcens_pois <- function(x, lambda) {
  1 - ppois(x - 1, lambda)
}

dcens_pois_upper <- function(x, lambda, upper) {
  ppois(upper, lambda) - ppois(x - 1, lambda)
}

x <- seq(1, 40)
mu <- 30
u <- 35

dc <- dcens_pois(x, 30)
dcu <- dcens_pois_upper(x, 30, upper = u)
d <- dpois(x, 30)

# dev.off()
# plot(x, dc, type = "h", ylim = c(0, 1))
# par(add = FALSE)
# plot(x, dc, type = "h", ylim = c(0, 1))

dat <- data.frame(
  x = c(x, x, x),
  dens = c(d, dc, dcu),
  type = c(
    rep("dpois", length(x)),
    rep("dcens_pois", length(x)),
    rep("dcens_pois_upper", length(x))
  )
)
library(ggplot2)

ggplot(dat, aes(x, dens))+ facet_wrap(~type) + geom_col() +
  geom_vline(xintercept = 15, col = "red", lwd = 2)

# ------------------------------------------------------------------

mu <- seq(1, 100)
x <- 30
u <- 60

dc <- dcens_pois(x, mu)
dcu <- dcens_pois_upper(x, mu, upper = u)
d <- dpois(x, mu)

dat <- data.frame(
  mu = c(mu, mu, mu),
  dens = c(d, dc, dcu),
  type = c(
    rep("dpois", length(mu)),
    rep("dcens_pois", length(mu)),
    rep("dcens_pois_upper", length(mu))
  )
)

g <- ggplot(dat, aes(mu, dens))+ facet_wrap(~type) + geom_col() +
  geom_vline(xintercept = x, col = "red", lwd = 2)

glog <- ggplot(dat, aes(mu, -log(dens)))+ facet_wrap(~type) + geom_line() +
  geom_vline(xintercept = x, col = "red", lwd = 2)
cowplot::plot_grid(g, glog, ncol = 1)


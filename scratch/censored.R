library(dplyr)
library(ggplot2)
theme_set(theme_light())

logspace_sub <- function(log_a, log_b) {
  if (log_a == log_b) {
    return(-Inf)  # Result should be zero in linear space, which is -Inf in log space
  } else if (log_a < log_b) {
    return(log_b + log1p(-exp(log_a - log_b)))
  } else {
    return(log_a + log1p(-exp(log_b - log_a)))
  }
}

dcens_pois <- function(x, lambda) {
  1 - ppois(x - 1, lambda)
}

dcens_pois_upper <- function(x, lambda, upper) {
  ppois(upper, lambda) - ppois(x - 1, lambda)
}

dcens_pois_joe <- function(x, lambda) {
  tmp_ll = log(ppois(x-1, lambda)) # // F(lower-1)
  tmp_ll = logspace_sub(0, tmp_ll) # // 1 - F(lower-1)
  exp(tmp_ll)
}

dcens_pois_upper_joe <- function(x, lambda, upper) {
  tmp_ll = log(ppois(upper, lambda)) # // F(upr)
  if (x > 0) {
    tmp_ll = logspace_sub(tmp_ll, log(ppois(x-1, lambda))) #// F(upr) - F(lwr-1) iff lwr>0
  }
  exp(tmp_ll)
}

dcens_pois_joe(9, 10)
dcens_pois(9, 10)

dcens_pois_upper(x = 9, lambda = 10, upper = 12)
dcens_pois_upper_joe(x = 9, lambda = 10, upper = 12)
dcens_pois_upper(x = 0, lambda = 10, upper = 12)
dcens_pois_upper_joe(x = 0, lambda = 10, upper = 12)

# -----------------------------

dcenspois_sean <- function(x, lambda, upr, give_log = 0L) {
  if (is.na(upr)) { # full right censored
    if (x == 0) {
      ll <- 0
    } else {
      ll <- log(dcens_pois_joe(x, lambda))
    }
  } else if (upr > x) { # upper truncated right censored
    ll <- log(dcens_pois_upper_joe(x, lambda, upr))
  } else if (x == upr) { # not censored
    ll <- dpois(x, lambda, log = TRUE)
  }
  if (give_log) {
    return(ll)
  } else {
    return(exp(ll))
  }
}

dcens_pois_upper(9, 10, 12)
dcenspois_sean(9, 10, upr = 12)

dcens_pois_upper(9, 10, 12)
dcenspois_sean(9, 10, upr = 12)

dcens_pois(9, 10)
dcenspois_sean(9, 10, upr = NA)

dcens_pois(0, 10)
dcenspois_sean(0, 10, upr = NA)

dcens_pois(0, 10)
dcenspois_sean(0, 10, upr = 11)

# ------------------



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


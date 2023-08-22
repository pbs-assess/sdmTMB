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
logspace_subtraction <- function(log_x, log_y) {
  max_log <- max(log_x, log_y)
  diff_exp <- exp(log_x - max_log) - exp(log_y - max_log)
  log_diff <- log(diff_exp) + max_log
  return(log_diff)
}
logspace_sub(-0.001, -0.00002)
logspace_subtraction(-0.001, -0.00002)

log_poisson_cdf <- function(x, lambda) {
  log_cdf <- sum(dpois(0:x, lambda, log = TRUE))
  return(log_cdf)
}
log_poisson_cdf(51, 1)
pp(51, 1)

x <- seq(1, 10)
plot(x, exp(ppois(x, 4, log.p = T)))

x <- seq(1, 10)
exp(log_poisson_cdf(1, 4))
ppois(1, 4)

#########


x1 <- dpois(0, 1)
x2=dpois(1, 1)
x3=dpois(2, 1)

log(x1 + x2 + x3)
ppois(2, 1, log.p = T)

x1 <- dpois(0, 1, log = T)
x2=dpois(1, 1, log = T)
x3=dpois(2, 1, log = T)
logspace_add(logspace_add(x1, x2), x3)

log(x1 + x2 + x3)
ppois(2, 1, log.p = T)




x1 <- dpois(0, 1)
x2=dpois(1, 1)
x3=dpois(2, 1)

x1 + x2 + x3
ppois(2, 1)



log_poisson_cdf <- function(q, lambda) {
  log_cdf <- 0.0
  for (k in 0:q) {
    log_cdf <- logspace_add(log_cdf, dpois(k, lambda, log = TRUE))
  }
  return(log_cdf)
}

# Helper function to add values in log space
logspace_add <- function(log_a, log_b) {
  if (log_a < log_b) {
    return(log_b + log1p(exp(log_a - log_b)))
  } else {
    return(log_a + log1p(exp(log_b - log_a)))
  }
}


# Test the function
q <- 2
lambda <- 1.5
log_poisson_cdf(q, lambda)
ppois(q, lambda)

#######



# double logPoissonCDF(int k, double lambda) {
#   double logSum = -lambda;
#   double logTerm = -lambda;
#
#   for (int i = 1; i <= k; i++) {
#     logTerm += log(lambda / i);
#     logSum = logSum + logTerm;
#   }
#
#   return exp(logSum);
# }

logPoissonCDF <- function(k, lambda) {
  logSum <- -lambda
  logTerm <- -lambda

  for (i in 1:k) {
    logTerm <- logTerm + log(lambda / i)
    logSum <- logSum + logTerm
  }

  return(logSum)
}

logPoissonCDF <- function(k, lambda) {
  logSum <- ppois(0, lambda, log.p = TRUE)
  for (i in 1:k) {
    logSum <- logSum + dpois(i, lambda, log = TRUE)
  }
  return(logSum)
}

logPoissonCDF <- function(k, lambda) {
  logSum <- 0
  logTerm <- 0

  for (i in 1:k) {
    logTerm <- log(lambda) + logTerm - log(i)
    logSum <- logSum + exp(logTerm)
  }

  return(log(logSum))
}


ppois(51, 1, log.p = TRUE)

logPoissonCDF(51, 1)

dcens_pois <- function(x, lambda) {
  1 - ppois(x - 1, lambda)
}

dcens_pois_upper <- function(x, lambda, upper) {
  ppois(upper, lambda) - ppois(x - 1, lambda)
}

dcens_pois_joe <- function(x, lambda) {
  tmp_ll = ppois(x-1, lambda, log.p = TRUE) # // F(lower-1)
  tmp_ll = logspace_sub(0, tmp_ll) # // 1 - F(lower-1)
  tmp_ll
}

dcens_pois_upper_joe <- function(x, lambda, upper) {
  tmp_ll = ppois(upper, lambda, log.p = TRUE) # // F(upr)
  if (x > 0) {
    tmp_ll = logspace_sub(tmp_ll, ppois(x-1, lambda, log.p = TRUE)) #// F(upr) - F(lwr-1) iff lwr>0
  }
  tmp_ll
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
      ll <- dcens_pois_joe(x, lambda)
    }
  } else if (upr > x) { # upper truncated right censored
    ll <- dcens_pois_upper_joe(x, lambda, upr)
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


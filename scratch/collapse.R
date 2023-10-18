library(RTMB)

set.seed(123)
x <- runif(200)
g <- rep(gl(5, 1), length(x))[seq_along(x)]
re <- rnorm(length(unique(g)), 0, sd = 0.01)
y <- rnorm(length(x), 0.5 + 0.3 * x + re[g], sd = 0.2)
plot(x, y)

dat <- list(x = x, y = y, g = as.integer(g), n = length(x))
par <- list(b0 = 0, b1 = 0, log_sigma = rep(0, 2L), theta = rep(0, length(re)))

jnll <- function(par) {
  getAll(par, dat)
  mu <- numeric(length = n)
  for (i in 1:n) {
    mu[i] <- b0 + x[i] * b1 + theta[g[i]]
  }
  nll <- 0
  nll <- nll - sum(dnorm(y, mu, exp(log_sigma[1]), log = TRUE))
  nll <- nll - sum(dnorm(theta, 0, exp(log_sigma[2]), log = TRUE))
  nll
}


# thresh <- log(1e-02)
# fixvals <- log(1e-03)
# re2zero <- "theta"

random <- "theta"
obj <- MakeADFun(jnll, par, random = random)

get_pars <- function(object) {
  ee <- object$env
  x <- ee$last.par.best
  if (length(ee$random) > 0) x <- x[-ee$random]
  p <- ee$parList(x = x)
  p
}

fit <- nlminb(obj$par, obj$fn, obj$gr)

tocheck <- list(list(log_sigma = c(0.01, 0.01)))
tomap <- list('theta')

pars <- get_pars(obj)

todo <- lapply(tocheck, function(.x) {
  nm <- names(.x)
  any(pars[[nm]] < .x[[nm]])
})

if (any(unlist(todo))) {
  m <- obj$env$map
  for (i in seq_along(todo)) {
    if (todo[[i]]) {
      nm <- names(tocheck[[i]])
      tc <- tocheck[[i]][[1]]
      tm <- tomap[[i]]

      # FIXME!!!
      m[[nm]] <- factor(c(1, NA))
      # pars[[names(tc)]] <- !is.na(tc[[1]])
      p <- pars[[tc]]
      pars[[tc]][p < tc] <- tc[p < tc]
      # FIXME!!! end

      lp <- length(pars[[tm]])
      m[[tm]] <- factor(rep(NA, lp))
      pars[[tm]] <- rep(0, lp)

      random <- random[!random %in% names(tm)]
    }
  }
  if (length(random) == 0L) random <- NULL
  obj <- MakeADFun(jnll, pars, map = m, random = random)
  fit <- nlminb(obj$par, obj$fn, obj$gr)
}

sdr <- sdreport(obj)
summary(sdr)

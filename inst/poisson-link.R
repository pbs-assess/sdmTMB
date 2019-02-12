# ------------------------------------------------------------
# Do it once:

set.seed(3854)
N <- 100
x_i <- rnorm(N)

b0_n <- -0.2
b1_n <- 0.5
n_i <- exp(b0_n + b1_n * x_i)

b0_w <- -0.1
b1_w <- 0.7
w_i <- exp(b0_w + b1_w * x_i)
a_i <- rep(0.4, N)

# need to get r_i and p_i
p_i <- 1 - exp(-a_i * n_i)
r_i <- (n_i / p_i) * w_i
p_i
r_i

plot(x_i, p_i)
plot(x_i, r_i)

# data:
y_present <- rbinom(N, 1, p_i)
y_present
CV <- 0.2
y_rate <- rgamma(N, shape = 1/(CV^2), rate = 1/((CV^2) * r_i))
y_rate <- y_rate * y_present

plot(x_i, y_rate)

nll_poisson_link <- function(par) {
  # a_i: area swept data
  # x_i: covariate data
  # y_rate: rate data (e.g. kg/hr)

  y_present <- ifelse(y_rate > 0, 1, 0)
  n_i <- exp(par[1] + par[2] * x_i)
  w_i <- exp(par[3] + par[4] * x_i)

  p_i <- 1 - exp(-a_i * n_i)
  r_i <- (n_i / p_i) * w_i

  nll1 <- -sum(dbinom(y_present, size = 1, prob = p_i, log = TRUE))
  nll2 <- -sum(dgamma(y_rate[as.logical(y_present)], shape = 1/(exp(par[5])^2),
    rate = 1/((exp(par[5])^2) * r_i[as.logical(y_present)]), log = TRUE))
  nll1 + nll2
}

a_i
x_i
y_present
y_rate
rm(y_present)

nll_poisson_link(c(b0_n, b1_n, b0_w, b1_w, CV))

start <- c(rep(0, 5))
m <- nlminb(start, nll_poisson_link)
m

# -------------------------------------------------
# Sim test;

set.seed(19254)
N <- 200
x_i <- rnorm(N)

b0_n <- -0.2
b1_n <- 0.5
n_i <- exp(b0_n + b1_n * x_i)

b0_w <- 0.1
b1_w <- 0.3
w_i <- exp(b0_w + b1_w * x_i)
a_i <- rep(1, N)

# need to get r_i and p_i
p_i <- 1 - exp(-a_i * n_i)
r_i <- (n_i / p_i) * w_i
p_i
r_i

plot(x_i, p_i)
plot(x_i, r_i)

CV <- 0.1

out <- lapply(seq_len(100), function(i) {
  set.seed(i)

  # data:
  y_present <- rbinom(N, 1, p_i)
  y_rate <- rgamma(N, shape = 1/(CV^2), rate = 1/((CV^2) * r_i))
  y_rate <- y_rate * y_present

  nll_poisson_link <- function(par) {
    # a_i: area swept data
    # x_i: covariate data
    # y_rate: rate data (e.g. kg/hr)
    n_i <- exp(par[1] + par[2] * x_i)
    w_i <- exp(par[3] + par[4] * x_i)

    p_i <- 1 - exp(-a_i * n_i)
    r_i <- (n_i / p_i) * w_i

    nll1 <- -sum(dbinom(y_present, size = 1, prob = p_i, log = TRUE))
    nll2 <- -sum(dgamma(y_rate[as.logical(y_present)], shape = 1/(exp(par[5])^2),
      rate = 1/((exp(par[5])^2) * r_i[as.logical(y_present)]), log = TRUE))
    nll1 + nll2
  }

  start <- c(rep(0, 5))
  m <- nlminb(start, nll_poisson_link)
  m$par
})

out <- plyr::ldply(out)
par(mfrow = c(2, 3))
hist(out[,1]);abline(v = b0_n, col = "red")
hist(out[,2]);abline(v = b1_n, col = "red")
hist(out[,3]);abline(v = b0_w, col = "red")
hist(out[,4]);abline(v = b1_w, col = "red")
hist(exp(out[,5]));abline(v = CV, col = "red")

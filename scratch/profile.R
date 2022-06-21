pl <- m1$tmb_obj$env$parList()

m1$tmb_obj$fn(m1$tmb_params)

m1$tmb_obj$env$parList()

# range = sqrt(8) / exp(ln_kappa)
# exp(ln_kappa) * range = sqrt(8)
# exp(ln_kappa) = sqrt(8) / range
# ln_kappa = log (sqrt(8) / range)

m4$tmb_obj$fn(m4$model$par)

ln_kappas <- seq(1, 40, length.out = 50)
# ln_kappas <- seq(1, 40, length.out = 50)
out <- sapply(ln_kappas, function(x) {
  .par <- m4$model$par
  .par[names(.par) == "ln_kappa"] <- log (sqrt(8) / x)
  as.numeric(m4$tmb_obj$fn(.par))
})

ln_kappas <- seq(-2, 0.3, length.out = 9)
ln_tau_Os <- seq(-0.5, 0.5, length.out = 10)

out <- sapply(ln_kappas, function(x) {
  sapply(ln_tau_Os, function(y) {
    .par <- m4$model$par
    .par[names(.par) == "ln_kappa"] <- x
    .par[names(.par) == "ln_tau_O"] <- y
    as.numeric(m4$tmb_obj$fn(.par))
  })})

# m4

# sigma_Os = 1 / sqrt(4 * 3.141 * exp(2 * ln_tau_Os + 2.0 * ln_kappa)

image(z = out, y = ln_kappas, x = ln_tau_Os)

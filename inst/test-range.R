x <- stats::runif(201, -1, 1)
y <- stats::runif(201, -1, 1)
kappa <- 2
sigma_O <- .3 # SD of spatial process
phi <- 0.1 # observation error

rf_sim <- function(model, x, y) {
  set.seed(sample.int(1e5L, 1L))
  set.seed(4)
  suppressMessages(
    RandomFields::RFsimulate(model = model, x = x, y = y)$variable1
  )
}

# .scale <- 1 / (sqrt(8) / kappa)
.scale <- 1/kappa
rf_omega <- RandomFields::RMmatern(nu = 1, var = sigma_O^2, scale = .scale)
omega_s <- rf_sim(model = rf_omega, x, y)
s <- data.frame(aa = rep(NA, length(x)))
s$x <- x
s$y <- y
s$z <- rnorm(length(omega_s), omega_s, phi)

library(ggplot2)
ggplot(s, aes(x, y, colour = z)) + geom_point() + scale_color_gradient2()

# spde <- sdmTMB::make_spde(x = s$x, y = s$y, n_knots = 100)
loc_xy <- cbind(s$x, s$y)
bnd <- INLA::inla.nonconvex.hull(as.matrix(loc_xy), convex = -0.3)
mesh <- INLA::inla.mesh.2d(
  boundary = bnd,
  max.edge = c(.1, 1),
  offset = -0.05,
  cutoff = c(.1, 1),
  min.angle = 1
)
spde <- sdmTMB::make_spde(s$x, s$y, mesh = mesh)
sdmTMB::plot_spde(spde)

tmb_data <- list(
  y_i = s$z,
  n_s = nrow(spde$mesh$loc),
  s_i = spde$cluster - 1L,
  spde = spde$spde$param.inla[c("M0", "M1", "M2")]
)

tmb_params <- list(
  ln_tau_O = 0,
  ln_kappa = 0,
  ln_phi = 0,
  omega_s = rep(0, length(omega_s))
)

library(TMB)
compile("inst/test_range.cpp")
dyn.load(dynlib("inst/test_range"))

tmb_random <- "omega_s"

tmb_obj1 <- TMB::MakeADFun(
  data = tmb_data, parameters = tmb_params,
  DLL = "test_range", random = tmb_random
)

tmb_opt1 <- stats::nlminb(
  start = tmb_obj1$par, objective = tmb_obj1$fn,
  gradient = tmb_obj1$gr
)

print(exp(tmb_opt1$par)[2])
r <- tmb_obj1$report()$sigma_O
print(r)

# sqrt(8) / exp(tmb_opt1$par)[[2]]

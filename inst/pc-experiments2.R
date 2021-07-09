# Experiments based on https://arxiv.org/abs/2103.09929

library(sdmTMB)
library(ggplot2)
d <- pcod
pcod_spde <- make_mesh(d, c("X", "Y"), cutoff = 20)
m <- sdmTMB(
  data = d, formula = density ~ 0 + as.factor(year) +
    depth_scaled + depth_scaled2,
  time = "year", spde = pcod_spde, family = tweedie(link = "log"),
  matern_prior_O = c(5, 0.05, 5, 0.05),
  matern_prior_E = c(5, 0.05, 5, 0.05),
  silent = FALSE
)

predict(m, newdata = qcs_grid)
SD0 <- TMB::sdreport(m$tmb_obj,
  getJointPrecision = TRUE
  # bias.correct = TRUE,
  # bias.correct.control = list(sd = TRUE)
)

mu <- c(SD0$par.fixed, SD0$par.random)

rmvnorm_prec <- function(mu, tmb_sd, n.sims) {
  z <- matrix(rnorm(length(mu) * n.sims), ncol = n.sims)
  L <- Matrix::Cholesky(tmb_sd[["jointPrecision"]], super = TRUE)
  z <- Matrix::solve(L, z, system = "Lt") ## z = Lt^-1 %*% z
  z <- Matrix::solve(L, z, system = "Pt") ## z = Pt    %*% z
  z <- as.matrix(z)
  mu + z
}

t.draws <- rmvnorm_prec(mu = mu, tmb_sd = SD0, n.sims = 300)
parnames <- c(names(SD0[["par.fixed"]]), names(SD0[["par.random"]]))
p <- list()
for (i in unique(parnames)) {
  p[[i]] <- t.draws[parnames == i, , drop = FALSE]
}

# ---------------------------------

lp <- m$tmb_obj$env$last.par.best

.lp <- c()
for (i in seq_along(p)) {

}
r <- m$tmb_obj$report(lp)

# ---------------------------------
object <- m
tmb_data <- object$tmb_data
newdata <- qcs_grid
xy_cols <- c("X", "Y")
newdata$sdm_orig_id <- seq(1, nrow(newdata))
fake_newdata <- unique(newdata[, xy_cols])
fake_newdata[["sdm_spatial_id"]] <- seq(1, nrow(fake_newdata)) - 1L
newdata <- base::merge(newdata, fake_newdata,
  by = xy_cols,
  all.x = TRUE, all.y = FALSE
)
newdata <- newdata[order(newdata$sdm_orig_id), , drop = FALSE]
proj_mesh <- INLA::inla.spde.make.A(object$spde$mesh,
  loc = as.matrix(fake_newdata[, xy_cols, drop = FALSE])
)

nd <- newdata
response <- sdmTMB:::get_response(object$formula)
# sdmTMB_fake_response <- FALSE
if (!response %in% names(nd)) {
  nd[[response]] <- 0 # fake for model.matrix
  # sdmTMB_fake_response <- TRUE
}
proj_X_ij <- model.matrix(object$formula, data = nd)

# dim(proj_mesh)
# dim(p$epsilon_st)
n_t <- m$tmb_data$n_t

eta_matrix <- matrix(nrow = nrow(newdata), ncol = ncol(t.draws))

for (j in seq_len(ncol(eta_matrix))) {
  # // Spatial and spatiotemporal random fields:
  proj_re_sp <- proj_mesh %*% p$omega_s[, j, drop = FALSE]
  proj_re_sp_st_all <- rep(proj_re_sp, n_t)
  proj_re_st_temp <- matrix(nrow = nrow(proj_mesh), ncol = n_t)
  proj_re_st <- matrix(nrow = nrow(proj_mesh), ncol = n_t)

  eps_st <- p$epsilon_st[, j, drop = TRUE]
  eps_st <- matrix(eps_st, ncol = n_t)
  for (i in 1:n_t) {
    proj_re_st_temp[, i] <- as.numeric(proj_mesh %*% eps_st[, i, drop = FALSE])
    proj_re_st[, i] <- proj_re_st_temp[, i]
  }

  proj_fe <- proj_X_ij %*% p$b_j[, j, drop = FALSE]

  # // Spatially varying coefficients:
  # proj_re_sp_trend <- rep(0, nrow(proj_X_ij))
  # proj_re_sp_slopes <- rep(0, nrow(proj_X_ij))

  # if (spatial_trend) {
  #   vector<Type> proj_re_sp_slopes_all = proj_mesh * omega_s_trend;
  #   for (int i = 0; i < proj_X_ij.rows(); i++) {
  #     proj_re_sp_trend(i) = proj_re_sp_slopes_all(proj_spatial_index(i)) * proj_t_i(i);
  #     proj_re_sp_slopes(i) = proj_re_sp_slopes_all(proj_spatial_index(i));
  #   }
  # }

  # // Pick out the appropriate spatial and/or or spatiotemporal values:
  # proj_re_st_vector <- rep(0, nrow(proj_X_ij))
  # proj_re_sp_st <- rep(0, nrow(proj_X_ij))
  # for (i 1:nrow(proj_X_ij)) {
  #   proj_re_sp_st(i) = proj_re_sp_st_all[proj_spatial_index[i]]
  #   proj_re_st_vector(i) = proj_re_st(proj_spatial_index[i], proj_year[i]]
  # }

  # just repeat it for now:

  proj_re_st_vector <- as.numeric(proj_re_st)
  proj_re_sp_st <- as.numeric(proj_re_sp)

  proj_rf <- proj_re_sp_st + proj_re_st_vector # + proj_re_sp_trend
  proj_eta <- proj_fe + proj_rf # // proj_fe includes RW and IID random effects

  eta_matrix[, j] <- proj_eta
}

eta <- apply(eta_matrix, 1, mean)
eta_med <- apply(eta_matrix, 1, median)
.sd <- apply(exp(eta_matrix), 1, sd)
cv <- .sd / exp(eta)
eta_se <- apply(eta_matrix, 1, sd)
lwr <- exp(eta - 2 * eta_se)
upr <- exp(eta + 2 * eta_se)

newdata$eta <- eta
newdata$eta_se <- eta_se
newdata$lwr <- lwr
newdata$upr <- upr
# dim(proj_eta)

# for comparison
p2 <- predict(m, newdata = qcs_grid)

# plot(p2$est - newdata$eta)
# plot(p2$est - eta_med)

ggplot(p2, aes(X, Y, fill = exp(est))) +
  geom_raster() +
  facet_wrap(~year) +
  scale_fill_viridis_c(option = "C", trans = "sqrt")

ggplot(newdata, aes(X, Y, fill = exp(eta))) +
  geom_raster() +
  facet_wrap(~year) +
  scale_fill_viridis_c(option = "C", trans = "sqrt")

ggplot(newdata, aes(X, Y, fill = cv)) +
  geom_raster() +
  facet_wrap(~year) +
  scale_fill_viridis_c(option = "B", trans = "log")

ggplot(newdata, aes(X, Y, fill = eta_se)) +
  geom_raster() +
  facet_wrap(~year) +
  scale_fill_viridis_c(option = "B", trans = "log")

ggplot(newdata, aes(X, Y, fill = lwr)) +
  geom_raster() +
  facet_wrap(~year) +
  scale_fill_viridis_c(option = "C", trans = "sqrt")

ggplot(newdata, aes(X, Y, fill = upr)) +
  geom_raster() +
  facet_wrap(~year) +
  scale_fill_viridis_c(option = "C", trans = "sqrt")

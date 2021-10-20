library(ggplot2)
library(tictoc)
library(sdmTMB)

pcod_spde <- make_mesh(pcod, c("X", "Y"), cutoff = 8)
plot(pcod_spde)
pcod_spde$mesh$n
m <- sdmTMB(density ~ 0 + as.factor(year),
  data = pcod, spde = pcod_spde, family = tweedie(link = "log"),
  time = "year", silent = FALSE)
print(m)


# m <- sdmTMB(density ~ 0 + as.factor(year),
#   data = pcod, spde = pcod_spde, family = tweedie(link = "log"),
#   time = "year", silent = FALSE)

# tic()
# o1 <- TMB::sdreport(m$tmb_obj)
# ot1 <- toc()
#
# tic()
# o2 <- TMB::sdreport(m$tmb_obj, getJointPrecision = TRUE)
# ot2 <- toc()
#
# ot1$toc - ot1$tic
# ot2$toc - ot2$tic
# pryr::object_size(o1)
# pryr::object_size(o2)

timings <- list()
tic()
p <- predict(m, newdata = qcs_grid, return_tmb_object = TRUE)
ind_tmb <- get_index(p)
timings[[1]] <- toc()

tic()
p <- predict(m, newdata = qcs_grid, return_tmb_object = TRUE)
ind_tmb_biascor <- get_index(p, bias_correct = TRUE)
timings[[2]] <- toc()

set.seed(1)
tic()
p <- predict(m, newdata = qcs_grid, sims = 400L)
# ind_sim <- get_index_sims(m, p, newdata = qcs_grid, est_function = function(x) exp(mean(log(x))))
ind_sim <- get_index_sims(p)
# ind_sim <- get_index_sims(p, return_sims = TRUE)
timings[[3]] <- toc()
# ind_sim <- mutate(ind_sim,
#   lwr = exp(log_est - 1.96 * se),
#   upr = exp(log_est + 1.96 * se))

.t <- lapply(timings, function(x) as.numeric(round(x$toc - x$tic, 1L)))

library(dplyr)
ind <- bind_rows(
  select(ind_tmb, year, est, lwr, upr) %>% mutate(type = paste0("TMB (", .t[[1]], " sec)")),
  select(ind_tmb_biascor, year, est, lwr, upr) %>% mutate(type = paste0("TMB bias corrected (", .t[[2]], " sec)")),
  select(ind_sim, year, est, lwr, upr) %>% mutate(type = paste0("Simulated (", .t[[3]], " sec)"))
)
ggplot(ind, aes(year, est, ymin = lwr, ymax = upr)) +
  geom_ribbon(alpha = 0.3, mapping = aes(fill = type)) +
  geom_line(aes(colour = type)) +
  scale_fill_brewer(palette = "Dark2") +
  scale_colour_brewer(palette = "Dark2")

pcod_spde <- make_mesh(pcod, c("X", "Y"), cutoff = 8)
plot(pcod_spde)
pcod_spde$mesh$n
m_tmb <- sdmTMB(density ~ 0 + as.factor(year),
  data = pcod, spde = pcod_spde, family = tweedie(link = "log"),
  time = "year", silent = FALSE, start = list(ln_kappa = -1.581, ln_tau_E = -0.146, ln_tau_O = -0.648), map = list(ln_kappa = factor(NA), ln_tau_E = factor(NA), ln_tau_O = factor(NA)))
print(m_tmb)
print(m)
m_tmb$tmb_obj$env$parList()$ln_tau_E
m_tmb$tmb_obj$env$parList()$ln_tau_O
m_tmb$tmb_obj$env$parList()$ln_kappa

library(tmbstan)
m_stan <- tmbstan(m_tmb$tmb_obj, silent = FALSE,
  iter = 200, chains = 1, control = list(adapt_delta = 0.9))
m_stan

length(m_tmb$tmb_obj$env$last.par.best)
names(m_tmb$tmb_obj$env$last.par.best)[1:100]
post <- rstan::extract(m_stan)
names(post)

# ---------------------------
# This is the important code to convert Stan posterior samples
# into a matrix where each column can be fed into `tmb_obj$report()`
post <- rstan::extract(m_stan)
p_names <- unique(names(m_tmb$tmb_obj$env$last.par.best))
post_matrix <- matrix(ncol = nrow(post$b_j),
  nrow = length(m_tmb$tmb_obj$env$last.par.best))
for (i in seq_len(ncol(post_matrix))) {
  post_pars <- list()
  for (j in seq_along(p_names)) {
    par_j <- p_names[j]
    if (is.matrix(post[[par_j]])) {
      post_pars[[j]] <- post[[par_j]][i, , drop = TRUE]
    } else {
      post_pars[[j]] <- post[[par_j]][i]
    }
  }
  post_matrix[, i] <- unlist(post_pars)
}
# ---------------------------

# post_matrix <- matrix(ncol = nrow(post$b_j),
#   nrow = length(m_tmb$tmb_obj$env$last.par.best))
# for (i in seq_len(ncol(post_matrix))) {
#   post_pars <- numeric()
#   post_pars <- c(post_pars, post$b_j[i,,drop=TRUE])
#   post_pars <- c(post_pars, post$thetaf[i])
#   post_pars <- c(post_pars, post$ln_phi[i])
#   post_pars <- c(post_pars, post$omega_s[i,,drop=TRUE])
#   post_pars <- c(post_pars, post$epsilon_st[i,,drop=TRUE])
#   stopifnot(length(post_pars) == length(m_tmb$tmb_obj$env$last.par.best))
#   post_matrix[,i] <- post_pars
# }

# post_tmb <- matrix(ncol = nrow(post$b_j), nrow = length(post_))
#   m_tmb$tmb_obj$report(post_pars)

# -------- prediction ------------------

# new: from above:
object <- m_tmb # new
tmb_data <- object$tmb_data
tmb_data$do_predict <- 1L

newdata <- qcs_grid
newdata$sdm_orig_id <- seq(1, nrow(newdata))
xy_cols <- c("X", "Y") # new
fake_newdata <- unique(newdata[,xy_cols])
fake_newdata[["sdm_spatial_id"]] <- seq(1, nrow(fake_newdata)) - 1L

newdata <- base::merge(newdata, fake_newdata, by = xy_cols,
  all.x = TRUE, all.y = FALSE)
newdata <- newdata[order(newdata$sdm_orig_id),, drop=FALSE]

proj_mesh <- INLA::inla.spde.make.A(object$spde$mesh,
  loc = as.matrix(fake_newdata[,xy_cols, drop = FALSE]))

# this formula has breakpt() etc. in it:
thresh <- check_and_parse_thresh_params(object$formula, newdata)
formula <- thresh$formula # this one does not

nd <- newdata
response <- get_response(object$formula)
sdmTMB_fake_response <- FALSE
if (!response %in% names(nd)) {
  nd[[response]] <- 0 # fake for model.matrix
  sdmTMB_fake_response <- TRUE
}

# deal with prediction IID random intercepts:
RE_names <- barnames(object$split_formula$reTrmFormulas)
## not checking so that not all factors need to be in prediction:
# fct_check <- vapply(RE_names, function(x) check_valid_factor_levels(data[[x]], .name = x), TRUE)
proj_RE_indexes <- vapply(RE_names, function(x) as.integer(nd[[x]]) - 1L, rep(1L, nrow(nd)))

if (!"mgcv" %in% names(object)) object[["mgcv"]] <- FALSE
proj_X_ij <- matrix(999)
if (!object$mgcv) {
  proj_X_ij <- tryCatch({model.matrix(object$formula, data = nd)},
    error = function(e) NA)
}
if (object$mgcv || identical(proj_X_ij, NA)) {
  proj_X_ij <- mgcv::predict.gam(object$mgcv_mod, type = "lpmatrix", newdata = nd)
}
if (!is.null(object$time_varying)) {
  proj_X_rw_ik <- model.matrix(object$time_varying, data = nd)
} else {
  proj_X_rw_ik <- matrix(0, ncol = 1, nrow = 1) # dummy
}

area <- 1 # new
if (length(area) != nrow(proj_X_ij) && length(area) != 1L) {
  stop("`area` should be of the same length as `nrow(newdata)` or of length 1.", call. = FALSE)
}

se_fit <- FALSE # new
pop_pred <- FALSE # new
exclude_RE <- FALSE # new
tmb_data$proj_X_threshold <- thresh$X_threshold
tmb_data$area_i <- if (length(area) == 1L) rep(area, nrow(proj_X_ij)) else area
tmb_data$proj_mesh <- proj_mesh
tmb_data$proj_X_ij <- proj_X_ij
tmb_data$proj_X_rw_ik <- proj_X_rw_ik
tmb_data$proj_RE_indexes <- proj_RE_indexes
tmb_data$proj_year <- make_year_i(nd[[object$time]])
tmb_data$proj_lon <- newdata[[xy_cols[[1]]]]
tmb_data$proj_lat <- newdata[[xy_cols[[2]]]]
tmb_data$calc_se <- as.integer(se_fit)
tmb_data$pop_pred <- as.integer(pop_pred)
tmb_data$exclude_RE <- exclude_RE
tmb_data$calc_time_totals <- as.integer(!se_fit)
tmb_data$proj_spatial_index <- newdata$sdm_spatial_id
tmb_data$proj_z_i <- as.numeric(newdata[[object$time]])
tmb_data$proj_z_i <- tmb_data$proj_t_i - mean(unique(tmb_data$proj_t_i)) # center on mean

epsilon_covariate <- rep(0, length(unique(newdata[[object$time]])))
if (tmb_data$est_epsilon_model) {
  # covariate vector dimensioned by number of time steps
  time_steps <- unique(newdata[[object$time]])
  for (i in seq_along(time_steps)) {
    epsilon_covariate[i] <- newdata[newdata[[object$time]] == time_steps[i],
      object$epsilon_predictor, drop = TRUE][[1]]
  }
}
tmb_data$epsilon_predictor <- epsilon_covariate

new_tmb_obj <- TMB::MakeADFun(
  data = tmb_data,
  parameters = object$tmb_obj$env$parList(),
  map = object$tmb_map,
  random = object$tmb_random,
  DLL = "sdmTMB",
  silent = TRUE
)

old_par <- object$model$par
# need to initialize the new TMB object once:
new_tmb_obj$fn(old_par)

sims <- 50L # new
if (sims > 0) {
  if (!"jointPrecision" %in% names(object$sd_report)) {
    message("Rerunning TMB::sdreport() with `getJointPrecision = TRUE`.")
    sd_report <- TMB::sdreport(object$tmb_obj, getJointPrecision = TRUE)
  } else {
    sd_report <- object$sd_report
  }
  t_draws <- rmvnorm_prec(mu = new_tmb_obj$env$last.par.best, # new: don't need this now
    tmb_sd = sd_report, n_sims = sims) # new: don't need this now
  r <- apply(t_draws, 2, new_tmb_obj$report) # new: don't need this now
  t_draws <- post_matrix # new
  r <- apply(t_draws, 2, new_tmb_obj$report)
  out <- lapply(r, `[[`, "proj_eta")
  out <- do.call("cbind", out)
  rownames(out) <- nd[[object$time]] # for use in index calcs
  attr(out, "time") <- object$time
  # return(out) # new: commented this out
}

# --------- end prediction

# --- now into index calcs

obj <- out # new, setup
.t <- as.numeric(rownames(obj))
yrs <- sort(unique(.t))
yr_indexes <- lapply(yrs, function(x) which(.t %in% x))
agg_function <- function(x) sum(exp(x)) # new
out1 <- lapply(yr_indexes, function(x) {
  apply(obj[x, , drop = FALSE], 2L, agg_function)
})
return_sims <- TRUE # new
if (return_sims) {
  out2 <- lapply(seq_along(out1), function(i) {
    ret <- data.frame(
      .time = yrs[i], .value = out1[[i]],
      .iteration = seq_along(out1[[i]])
    )
    stats::setNames(ret, c(attr(obj, "time"), ".value", ".iteration"))
  })
  # return(do.call("rbind", out2))
  out2 <- do.call("rbind", out2)

  p <- predict(m_tmb, newdata = qcs_grid, sims = 200)
  index <- get_index_sims(p, return_sims = TRUE)

  bind_rows(mutate(out2, type = "Stan"), mutate(index, type = "TMB simulated")) %>%
    ggplot(aes(as.factor(year), .value, fill = type)) + geom_violin()


  ##
  ## x <- get_index_sims(m, p)
  ## ggplot(x, aes(year, est, ymin = lwr, ymax = upr)) + geom_line() + geom_ribbon(alpha = 0.4)
  ## x <- get_index_sims(m, p, return_sims = TRUE)
  ## ggplot(x, aes(as.factor(year), .value)) + geom_violin()

  # ---------------------------------------
  ## p_nobias <- predict(m, newdata = qcs_grid, sim = 2000)
  ## dim(p)
  ## #
  ## p2 <- predict(m, newdata = qcs_grid, return_tmb_object = TRUE)
  ## .i <- get_index(p2)
  ## .i_bias <- get_index(p2, bias_correct = TRUE)
  ##
  ##
  ## .m1 <- apply(p, 1, mean)
  ## plot(.m1, p2$data$est);abline(a = 0, b = 1)
  ##
  ##
  ## .mean <- numeric(length = length(unique(qcs_grid$year)))
  ## .lwr <- numeric(length = length(unique(qcs_grid$year)))
  ## .upr <- numeric(length = length(unique(qcs_grid$year)))
  ## for (i in seq_along(unique(qcs_grid$year))) {
  ##   .year <- unique(qcs_grid$year)[i]
  ##   y2003 <- p[qcs_grid$year == .year, , drop=FALSE]
  ##   y2003_sum <- apply(y2003, 2, function(x) sum(exp(x)))
  ##   .mean[i] <- mean(y2003_sum)
  ##   .lwr[i] <- quantile(y2003_sum, probs = 0.025)
  ##   .upr[i] <- quantile(y2003_sum, probs = 0.975)
  ## }
  ##
  ## .mean <- numeric(length = length(unique(qcs_grid$year)))
  ## .lwr <- numeric(length = length(unique(qcs_grid$year)))
  ## .upr <- numeric(length = length(unique(qcs_grid$year)))
  ## for (i in seq_along(unique(qcs_grid$year))) {
  ##   .year <- unique(qcs_grid$year)[i]
  ##   y2003 <- p[qcs_grid$year == .year, , drop=FALSE]
  ##   y2003_sum <- apply(y2003, 2, function(x) sum(exp(x)))
  ##   .mean[i] <- mean(y2003_sum)
  ##   .lwr[i] <- quantile(y2003_sum, probs = 0.025)
  ##   .upr[i] <- quantile(y2003_sum, probs = 0.975)
  ## }
  ##
  ## .mean_nobias <- numeric(length = length(unique(qcs_grid$year)))
  ## .lwr_nobias <- numeric(length = length(unique(qcs_grid$year)))
  ## .upr_nobias <- numeric(length = length(unique(qcs_grid$year)))
  ## for (i in seq_along(unique(qcs_grid$year))) {
  ##   .year <- unique(qcs_grid$year)[i]
  ##   y2003 <- p_nobias[qcs_grid$year == .year, , drop=FALSE]
  ##   y2003_sum <- apply(y2003, 2, function(x) sum(exp(x)))
  ##   .mean_nobias[i] <- mean(y2003_sum)
  ##   .lwr_nobias[i] <- quantile(y2003_sum, probs = 0.025)
  ##   .upr_nobias[i] <- quantile(y2003_sum, probs = 0.975)
  ## }
  ##
  ## plot(.mean, .i$est, asp = 1);abline(a = 0, b=1)
  ## plot(log(.mean), log(.i$est), asp = 1);abline(a = 0, b = 1)
  ## .mean
  ## .i$est
  ##
  ## cols <- RColorBrewer::brewer.pal(4, "Dark2")
  ## plot(.i$year, .i$est, type = "l", ylim = range(.i$lwr, .i$upr, .lwr_nobias, .lwr, .upr_nobias, .upr), col = cols[1], xlab = "Year", ylab = "Biomass density index")
  ## lines(.i$year, .i$lwr, type = "l", lty = 2, col = cols[1])
  ## lines(.i$year, .i$upr, type = "l", lty = 2, col = cols[1])
  ##
  ## lines(.i$year, .mean, col = cols[2])
  ## lines(.i$year, .lwr, col = cols[2], lty = 2)
  ## lines(.i$year, .upr, col = cols[2], lty = 2)
  ##
  ## lines(.i$year, .mean_nobias, col = cols[3])
  ## lines(.i$year, .lwr_nobias, col = cols[3], lty = 2)
  ## lines(.i$year, .upr_nobias, col = cols[3], lty = 2)
  ##
  ## lines(.i_bias$year, .i_bias$est, col = cols[4])
  ## lines(.i_bias$year, .i_bias$lwr, col = cols[4], lty = 2)
  ## lines(.i_bias$year, .i_bias$upr, col = cols[4], lty = 2)
  ##
  ## legend(2006, 5e5, col = cols, legend = c("TMB sum", "Sim", "Sim w bias correction", "TMB sum with bias correction"), lty = 1, bg = NA, bty = "n")
  ##
  ## # .m <- apply(y2003, 1, mean)
  ## # q1 <- apply(y2003, 1, quantile, probs = 0.025)
  ## # q2 <- apply(y2003, 1, quantile, probs = 0.975)
  ## hist(log(y2003_sum))
  ## abline(v = .i$log_est[.i$year == .year])
  ##
  ## plot(.i$lwr, q1)
  ##
  ##

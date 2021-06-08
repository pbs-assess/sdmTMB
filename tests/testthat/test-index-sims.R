## library(ggplot2)
## library(tictoc)
## library(sdmTMB)
##
## pcod_spde <- make_mesh(pcod, c("X", "Y"), cutoff = 30)
## plot(pcod_spde)
## pcod_spde$mesh$n
## m <- sdmTMB(density ~ 0 + as.factor(year),
##   data = pcod, spde = pcod_spde, family = tweedie(link = "log"),
##   time = "year")
## print(m)
##
## timings <- list()
## tic()
## p <- predict(m, newdata = qcs_grid, return_tmb_object = TRUE)
## ind_tmb <- get_index(p)
## timings[[1]] <- toc()
##
## tic()
## p <- predict(m, newdata = qcs_grid, return_tmb_object = TRUE)
## ind_tmb_biascor <- get_index(p, bias_correct = TRUE)
## timings[[2]] <- toc()
##
## tic()
## p <- predict(m, newdata = qcs_grid, sims = 100L)
## # ind_sim <- get_index_sims(m, p, newdata = qcs_grid, est_function = function(x) exp(mean(log(x))))
## ind_sim <- get_index_sims(p)
## ind_sim <- get_index_sims(p, return_sims = TRUE)
## timings[[3]] <- toc()
## # ind_sim <- mutate(ind_sim,
## #   lwr = exp(log_est - 1.96 * se),
## #   upr = exp(log_est + 1.96 * se))

## .t <- lapply(timings, function(x) as.numeric(round(x$toc - x$tic, 1L)))
##
##
## ind <- bind_rows(
##   select(ind_tmb, year, est, lwr, upr) %>% mutate(type = paste0("TMB (", .t[[1]], " sec)")),
##   select(ind_tmb_biascor, year, est, lwr, upr) %>% mutate(type = paste0("TMB bias corrected (", .t[[2]], " sec)")),
##   select(ind_sim, year, est, lwr, upr) %>% mutate(type = paste0("Simulated (", .t[[3]], " sec)"))
## )
## ggplot(ind, aes(year, est, ymin = lwr, ymax = upr)) +
##   geom_ribbon(alpha = 0.3, mapping = aes(fill = type)) +
##   geom_line(aes(colour = type)) +
##   scale_fill_brewer(palette = "Dark2") +
##   scale_colour_brewer(palette = "Dark2")
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

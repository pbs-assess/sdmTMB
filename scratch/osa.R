pcod_gaus <- subset(pcod, density > 0 & year >= 2013)
pcod_spde_gaus <- make_mesh(pcod_gaus, c("X", "Y"), cutoff = 15)
m_pos <- sdmTMB(log(density) ~ 0 + as.factor(year), time = "year",
  data = pcod_gaus, mesh = pcod_spde_gaus)

m_pos$tmb_data$sim_re

s <- m_pos$tmb_obj$simulate()
names(s)

r <- as.list(m_pos$sd_report, "Estimate")
plot(s$epsilon_st_A_vec)

plot(s$epsilon_st[,1], r$epsilon_st[,1]);abline(a = 0, b = 1)

plot(s$y_i, r$mu);abline(a = 0, b = 1)

plot(s$y_i, log(pcod_gaus$density))

plot(s$epsilon_st_A_vec, r$epsilon_st)

pcod_gaus$sim <- s$y_i

library(ggplot2)
ggplot(pcod_gaus, aes(X, Y, size = sim, colour = sim)) + geom_point() +
  scale_size_area() + scale_colour_viridis_c()

ggplot(pcod_gaus, aes(X, Y, size = log(density), colour = log(density))) + geom_point() +
  scale_size_area() + scale_colour_viridis_c()


r <- m_pos$tmb_obj$report()
pcod_gaus[["omega"]] <- r$omega_s_A
plot(r$omega_s_A, s$omega_s_A);abline(a = 0, b = 1)
m_pos$tmb_data$sim_re

s <- m_pos$tmb_obj$simulate(complete = TRUE)

head(s)
hist(s$y_i)
hist(log(pcod_gaus$density))


print(m_pos)
r <- as.list(m_pos$sd_report, "Estimate")
r$omega_s[1:10]

m_pos2 <- sdmTMB(log(density) ~ 0 + as.factor(year),
  data = pcod_gaus, mesh = pcod_spde_gaus)
print(m_pos2)
r2 <- as.list(m_pos2$sd_report, "Estimate")
r2$omega_s[1:10]

plot(r$omega_s, r2$omega_s);abline(a = 0, b =1)

m_pos3 <- sdmTMB(log(density) ~ 0 + as.factor(year),
  data = pcod_gaus, mesh = pcod_spde_gaus)
print(m_pos3)
r3 <- as.list(m_pos2$sd_report, "Estimate")
r2$omega_s[1:10]
plot(r3$omega_s, r2$omega_s);abline(a = 0, b =1)


o <- TMB::oneStepPredict(m_pos$tmb_obj, data.term.indicator = "keep",
  observation.name = "y_i", method = "oneStepGaussian")

pcod_spde <- make_mesh(pcod, c("X", "Y"), cutoff = 15)
m2 <- sdmTMB(density ~ 0 + as.factor(year),
  data = pcod, mesh = pcod_spde, family = tweedie())
print(m2)

o <- TMB::oneStepPredict(m2$tmb_obj, data.term.indicator = "keep",
  observation.name = "y_i", method = "oneStepGaussian", discrete = FALSE,
  range = c(0, max(pcod$density) * 1.1))

head(o)

r <- residuals(m_pos)
par(mfrow = c(1, 2))
qqnorm(r);qqline(r);abline(a = 0, b = 1, lty = 2)
qqnorm(o$residual);qqline(o$residual);abline(a = 0, b = 1, lty = 2)


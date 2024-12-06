library(lme4)
library(glmmTMB)

m <- glmmTMB(Reaction ~ Days + (Days || Subject), sleepstudy)
m2 <- sdmTMB(
  Reaction ~ Days + (1 | Subject), data = sleepstudy,
  experimental = list(slope_group = "Subject", slope_covariate = "Days"),
  spatial = "off")

m2$sd_report

logLik(m)
logLik(m2)

rs_glmmTMB <- ranef(m)[[1]][[1]][,2]
pars <- as.list(m2$sd_report, "Estimate")
rs_sdmTMB <- pars$RS

plot(rs_glmmTMB, rs_sdmTMB);abline(0, 1)

rs_glmmTMB - rs_sdmTMB

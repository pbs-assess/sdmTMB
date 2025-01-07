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

m2$tmb_data$RS_indexes
m2$tmb_data$RS_x

m2$tmb_data$proj_RS_indexes
m2$tmb_data$proj_RS_x
m2$slope_group
m2$slope_covariate

nd <- dplyr::select(sleepstudy, Subject) |> dplyr::distinct()
nd$Days <- 8
p <- predict(m, newdata = nd)
p2 <- predict(m2, newdata = nd)


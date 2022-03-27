d <- pcod_2011
mesh <- make_mesh(d, c("X", "Y"), cutoff = 20)

m <- sdmTMB(
  # data = d, formula = density ~ 1 + poly(log(depth), 2),
  data = d, formula = density ~ 1 + depth,
  mesh = mesh, family = sdmTMB::tweedie(link = "log"), spatial = "off"
)
m
visreg::visreg(m, "depth", scale = "response")

nd <- pcod_2011[1:3,,drop=FALSE]
p <- predict(m, newdata = nd)
p <- p$est

load_all("../visreg/")

family.sdmTMB <- function (object, ...) {
  object$family
}

# visreg::visreg(m, "depth", gg = TRUE)

library(glmmTMB)

m2 <- glmmTMB(
  # data = d, formula = density ~ 1 + poly(log(depth), 2),
  data = d, formula = density ~ 1 + depth,
  family = glmmTMB::tweedie(link = "log")
)
summary(m2)
visreg::visreg(m2, "depth", scale = "response")

p2 <- predict(m2, newdata = nd)
p
p2

pcod_gaus <- subset(pcod_2011, density > 0 & year >= 2013)
m_pos <- sdmTMB(density ~ 1 + depth_scaled + depth_scaled2,
  data = pcod_gaus, spatial = "off", family = Gamma(link = "log"))
m_pos2 <- glmmTMB(density ~ 1 + depth_scaled + depth_scaled2,
  data = pcod_gaus,family = Gamma(link = "log"))

nd$sdmTMB_X_ <- 500
nd$sdmTMB_Y_ <- 500
p <- predict(m_pos, newdata = nd)
p2 <- predict(m_pos2, newdata = nd)

p
p2


m_pos <- sdmTMB(density ~ 1 + depth_scaled + depth_scaled2,
  data = pcod, spatial = "off", family = tweedie(link = "log"))
m_pos2 <- glmmTMB(density ~ 1 + depth_scaled + depth_scaled2,
  data = pcod, family = glmmTMB::tweedie(link = "log"))

p <- predict(m_pos, newdata = nd)
p2 <- predict(m_pos2, newdata = nd)

p
p2

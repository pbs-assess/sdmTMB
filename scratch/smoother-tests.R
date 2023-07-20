library("sdmTMB")
library("lme4")
# library("glmmTMB")
library("gamm4")

# data("sleepstudy", package = "lme4")
# m1a <- lmer(Reaction ~ Days + (1 | Subject), data = sleepstudy, REML = TRUE)
# m1b <- sdmTMB(Reaction ~ Days + (1 | Subject), data = sleepstudy, spatial = "off", reml = TRUE)
#
# logLik(m1a)
# logLik(m1b)
#
# data(mcycle, package = 'MASS')
#
# detach("package:glmmTMB", unload=TRUE)
# m2b <- sdmTMB(accel ~ s(times), data = mcycle, spatial = "off", reml = FALSE)
# m2a <- gamm4(accel ~ s(times), data = mcycle, REML = FALSE)
# load_all("../glmmTMB/glmmTMB/")
# m2c <- glmmTMB(accel ~ s(times), data = mcycle, REML = FALSE)

data("sleepstudy", package = "lme4")
gamm1 <- gamm4(Reaction ~ s(Days), data = sleepstudy, REML = FALSE)
load_all("../glmmTMB/glmmTMB/")
gtmb1 <- glmmTMB(Reaction ~ s(Days), data = sleepstudy, REML = FALSE)
detach("package:glmmTMB", unload=TRUE)
load_all(".")
stmb1 <- sdmTMB(Reaction ~ s(Days), data = sleepstudy, reml = FALSE, spatial = "off")

logLik(gamm1$mer)
logLik(gtmb1)
logLik(stmb1)

coef(gamm1$mer)

gamm1$mer@beta
stmb1$model$par

gamm1$mer@frame

head(gtmb1$modelInfo$reTrms$cond$smooth_info[[1]]$re$Xf)
head(d$Xs[,1])

head(d$Zs_1_1)

######
plot(d$Xs[,1], gtmb1$modelInfo$reTrms$cond$smooth_info[[1]]$re$Xf[,1])
abline(0, 12.4)
######




est <- as.list(stmb1$sd_report, "Estimate", report = FALSE)
est$b_smooth
est$b_j
est$bs


library(brms)
brm1 <- brm(Reaction ~ s(Days), data = sleepstudy, chains = 1, iter = 50000, backend = "cmdstanr")
brm1

d <- brms::standata(brm1)
expect_equal(d$Xs[,1], stmb1$tmb_data$Xs[,1])
expect_equal(d$Zs_1_1, stmb1$tmb_data$Zs[[1]])
stmb1$model$par
exp(stmb1$model$par[["ln_smooth_sigma"]])

data(mcycle, package = 'MASS')
gamm2 <- gamm4(accel ~ s(times), data = mcycle, REML = FALSE)
load_all("../glmmTMB/glmmTMB/")
gtmb2 <- glmmTMB(accel ~ s(times), data = mcycle, REML = FALSE)
detach("package:glmmTMB", unload=TRUE)
load_all(".")
stmb2 <- sdmTMB(accel ~ s(times), data = mcycle, reml = FALSE, spatial = "off")

logLik(gamm2$mer)
gtmb2$fit$objective
logLik(stmb2)
fixef(gtmb2)$cond

stmb2$model$par
stmb2$tmb_random
est <- as.list(stmb2$sd_report, "Estimate", report = FALSE)
est$b_smooth
est$b_j
est$bs


#####################################

gamm3 <- gamm4(Reaction ~ s(Days), data = sleepstudy, REML = FALSE)
load_all("../glmmTMB/glmmTMB/")
gtmb3 <- glmmTMB(Reaction ~ s(Days), data = sleepstudy, REML = FALSE)
detach("package:glmmTMB", unload=TRUE)
load_all(".")
stmb3 <- sdmTMB(Reaction ~ s(Days), data = sleepstudy, reml = FALSE, spatial = "off")
stmb3

logLik(gamm3$mer)
gtmb3$fit$objective
logLik(stmb3)

fixef(gtmb3)$cond

stmb3$model$par
stmb3$tmb_random
est <- as.list(stmb3$sd_report, "Estimate", report = FALSE)
est$b_smooth
est$b_j
est$bs

########################

logLik(m2a$mer)
logLik(m2b)
logLik(m2c)
m2c$fit$objective

m2a$gam$smooth

m2a$mer@beta
m2b

m3a <- gamm4(accel ~ s(times), data = mcycle, REML = TRUE)
m3b <- sdmTMB(accel ~ s(times), data = mcycle, spatial = "off", reml = TRUE)

load_all("../glmmTMB/glmmTMB/")
m3c <- glmmTMB(accel ~ s(times), data = mcycle, REML = TRUE)


m3c$obj$env$random


m3b$tmb_random

logLik(m3a$mer)
logLik(m3b)
m3c$fit$objective

p2a <- predict(m2a$mer)
p2a
p2b <- predict(m2b)

plot(p2a, p2b$est)
cor(p2a, p2b$est)

p3a <- predict(m3a$mer)
p3a
p3b <- predict(m3b)

plot(p3a, p3b$est)
cor(p3a, p3b$est)

cor(p3a, p2a)
plot(p2a, p3a)
#


devtools::load_all("~/src/glmmTMB/glmmTMB/")

data("sleepstudy", package = "lme4")
m1 <- glmmTMB(formula = Reaction ~ s(Days), data = sleepstudy, REML = TRUE)
m1a <- glmmTMB(formula = Reaction ~ s(Days), data = sleepstudy, REML = TRUE)
logLik(m1)
logLik(m1a)


devtools::load_all("~/src/glmmTMB/glmmTMB/")
m2 <- glmmTMB(Reaction ~ s(Days, m = 2), data = sleepstudy, REML = FALSE)
m2a <- mgcv::gam(Reaction ~ s(Days, m = 2), data = sleepstudy, REML = FALSE)
detach("package:glmmTMB", unload=TRUE)
m2b <- sdmTMB(Reaction ~ s(Days, m = 2), data = sleepstudy, reml = FALSE, spatial = 'off')

logLik(m2)
logLik(m2a)
logLik(m2b)


library(mgcv)
library(gamm4)
set.seed(2)
dat <- gamSim(1,n=400,dist="normal",scale=2)
dat$f4 <- gl(15, 30, length = nrow(dat))

detach("package:glmmTMB", unload=TRUE)

b <- gamm4(y~s(x0) + (1|f4),data=dat, REML = FALSE)
devtools::load_all("~/src/glmmTMB/glmmTMB/")
b1 <- glmmTMB(y~s(x0),data=dat, REML = FALSE)
logLik(b)
logLik(b1)


devtools::load_all("~/src/glmmTMB/glmmTMB/")
m2 <- glmmTMB(Reaction ~ s(Days, m = 3), data = sleepstudy, REML = FALSE)
m2a <- mgcv::gam(Reaction ~ s(Days, m = 3), data = sleepstudy, REML = FALSE)
detach("package:glmmTMB", unload=TRUE)
m2b <- sdmTMB(Reaction ~ s(Days, m = 3), data = sleepstudy, reml = FALSE, spatial = 'off')

logLik(m2)
logLik(m2a)
logLik(m2b)


m2a <- gamm4::gamm4(Reaction ~ s(Days) + (1 | Subject), data = sleepstudy, REML = TRUE)


devtools::load_all("~/src/glmmTMB/glmmTMB/")
library(gamm4)
set.seed(0)
dat <- mgcv::gamSim(1, n = 400, scale = 2)
dat$fac <- fac <- as.factor(sample(1:20, 400, replace = TRUE))
dat$y <- dat$y + model.matrix(~ fac - 1) %*% rnorm(20) * .5
dat$y <- as.numeric(dat$y)

## defaults should be REML = FALSE
## doesn't match!?
br <- gamm4(y ~ s(x1), data = dat, random = ~ (1 | fac))
br1 <- glmmTMB(y ~ s(x0) + (1 | fac), data = dat)
logLik(br$mer) # looks like REML = TRUE
logLik(br1)

## but this matches:
br <- gamm4(y ~ s(x0), data = dat, random = ~ (1 | fac), REML = FALSE)
br1 <- glmmTMB(y ~ s(x0) + (1 | fac), data = dat, REML = FALSE)
logLik(br$mer)
logLik(br1)

## and so does this:
br <- gamm4(y ~ s(x0), data = dat, random = ~ (1 | fac), REML = TRUE)
br1 <- glmmTMB(y ~ s(x0) + (1 | fac), data = dat, REML = TRUE)
logLik(br$mer)
logLik(br1)

###############


br <- gamm4(y ~ s(x0, m = 3), data = dat, random = ~ (1 | fac), REML = TRUE)
devtools::load_all("~/src/glmmTMB/glmmTMB/")
br1 <- glmmTMB(y ~ s(x0, m = 3) + (1 | fac), data = dat, REML = TRUE)
detach("package:glmmTMB", unload=TRUE)
load_all(".")
br2 <- sdmTMB(y ~ s(x0, m = 3) + (1 | fac), data = dat, reml = TRUE, spatial = "off")

logLik(br$mer)
logLik(br1)
logLik(br2)

p0 <- predict(br$mer)
p1 <- predict(br1)
p2 <- predict(br2)
plot(p0, p1)
plot(p0, p2$est)

###################

br <- gam(y ~ s(x0), data = dat, REML = FALSE, method = "ML")
devtools::load_all("~/src/glmmTMB/glmmTMB/")
br1 <- glmmTMB(y ~ s(x0), data = dat, REML = FALSE)
detach("package:glmmTMB", unload=TRUE)
load_all(".")
br2 <- sdmTMB(y ~ s(x0), data = dat, reml = FALSE, spatial = "off")

logLik(br)
logLik(br1)
logLik(br2)

p0 <- predict(br)
p1 <- predict(br1)
p2 <- predict(br2)
plot(p0, p1)
plot(p0, p2$est)

##############

br <- gam(y ~ s(x0, m = 4), data = dat, REML = FALSE, method = "ML")
br <- gam(y ~ s(x0, m = 4), data = dat)
devtools::load_all("~/src/glmmTMB/glmmTMB/")
br1 <- glmmTMB(y ~ s(x0, m = 4), data = dat, REML = FALSE)
detach("package:glmmTMB", unload=TRUE)
detach("package:sdmTMB", unload=TRUE)
load_all(".")
br2 <- sdmTMB(y ~ s(x0, m = 4), data = dat, reml = FALSE, spatial = "off")

logLik(br)
logLik(br1)
logLik(br2)

p0 <- predict(br)
p1 <- predict(br1)
p2 <- predict(br2)
plot(p0, p1)
plot(p0, p2$est)

##############


detach("package:sdmTMB", unload=TRUE)
load_all("../glmmTMB/glmmTMB/")
m_glmmTMB <- glmmTMB(density ~ s(depth_scaled, m = 2) + (1|fy), data = d, REML = FALSE)



br <- gam(y ~ s(x0), data = dat, random = ~ (1 | fac), REML = TRUE)
br1 <- glmmTMB(y ~ s(x0) + (1 | fac), data = dat, REML = TRUE)
logLik(br$mer)
logLik(br1)

br <- gam(y ~ s(x0), data = dat)
library(sdmTMB)
m <- gam(y ~ s(x0), data = dat, REML = FALSE, method = "ML")


br <- gam(y ~ s(x0), data = dat, REML = FALSE, method = "ML")
br1 <- glmmTMB(y ~ s(x0), data = dat, REML = FALSE)
logLik(br)
logLik(br1)

br <- gam(y ~ s(x0), data = dat, REML = TRUE)
br1 <- glmmTMB(y ~ s(x0), data = dat, REML = TRUE)
logLik(br)
logLik(br1)
















detach("package:glmmTMB", unload = TRUE)
br2 <- sdmTMB(y ~ s(x0, m = 2) + (1 | fac), data = dat, reml = TRUE, spatial = "off")
logLik(br2)

b3 <- gam(y~s(x0, m = 2) + (1|fac),data=dat, reml = TRUE)

p0 <- predict(br$mer)
p1 <- predict(br1)
p2 <- predict(br2)
plot(p1, p0)
plot(p1, p2$est)






br <- gamm4(y~s(x0, m = 4),data=dat,random=~(1|fac), REML = TRUE)
devtools::load_all("~/src/glmmTMB/glmmTMB/")
br1 <- glmmTMB(y~s(x0, m = 4) + (1|fac),data=dat, REML = TRUE)
logLik(br$mer)
logLik(br1)

br <- gam(y~s(x0),data=dat, REML = TRUE)
br1 <- glmmTMB(y~s(x0),data=dat, REML = TRUE)
logLik(br)
logLik(br1)



br <- gam(y~s(x0, m = 4),data=dat, REML = TRUE)
br1 <- glmmTMB(y~s(x0, m = 4),data=dat, REML = TRUE)
logLik(br)
logLik(br1)





logLik(m2)
logLik(m2a)

m3 <- glmmTMB(Reaction ~ s(Days, bs = "re"), data = sleepstudy, REML = TRUE)

m2 <- gamm4::gamm4(Reaction ~ s(Days, fx = TRUE), data = sleepstudy, REML = TRUE)


sleepstudy$x <- rnorm(nrow(sleepstudy))
m2 <- glmmTMB(Reaction ~ s(Days, x), data = sleepstudy, REML = TRUE)
##########


library(glmmTMB)
devtools::load_all("~/src/glmmTMB/glmmTMB/")
data("sleepstudy", package = "lme4")

## matches
m1 <- glmmTMB(Reaction ~ s(Days), data = sleepstudy, REML = TRUE)
m1a <- mgcv::gam(Reaction ~ s(Days), data = sleepstudy, REML = TRUE, method = "ML")
logLik(m1)
#> 'log Lik.' -945.7767 (df=4)
logLik(m1a)
#> 'log Lik.' -945.7767 (df=4)

## setting m causes a mismatch
m2 <- glmmTMB(Reaction ~ s(Days, m = 2), data = sleepstudy, REML = TRUE)
m2a <- mgcv::gam(Reaction ~ s(Days, m = 2), data = sleepstudy, REML = TRUE)
logLik(m2)
#> 'log Lik.' -945.7767 (df=4)
logLik(m2a)
#> 'log Lik.' -950.1465 (df=3)

## arguments that are going to fail with unhelpful messages
m3 <- glmmTMB(Reaction ~ s(Days, bs = "re"), data = sleepstudy, REML = TRUE)
#> Error in dimnames(x) <- dn: length of 'dimnames' [2] not equal to array extent

m4 <- glmmTMB(Reaction ~ s(Days, fx = TRUE), data = sleepstudy, REML = TRUE)
#> Error in t.default(s$re$rand$Xr): argument is not a matrix

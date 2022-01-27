library(TMB)
library(MASS)
library(tweedie)
library(sdmTMB)

### Parameterization #######
# specify true values of metabolic index
Eo <- 0.4
n <- -0.29
Ao <- 200
# specify tweedie parameters
phi <- 5
s1 <- 1.5

SEED <- 1234
set.seed(SEED)

# set variance / covariance matrix of logpO2 and logtemp
sigma <- matrix( c(0.9204363, 0.23806646,
                   0.2380665, 0.09876192),
                 nrow = 2,
                 ncol = 2,
                 byrow = F)

# log means of po2 and temp
logmu <- c(-3.627865,  1.878259)
ndata <- 100
# specify true values of threshold and slope
b_threshold <- c(1.5, 4)
M = 1000 # average body mass (g)

#### END OF PARAMETER SPECIFICATION ###

### Create Fake Data ####
env.data  <- exp(mvrnorm(n = ndata, logmu,sigma))
# create dataframe
thedata <- data.frame(po2 = env.data[,1],
                      temp = env.data[,2],
                      mi = NA,
                      y = NA)


linear_threshold <- function(x, b_threshold) {
  if (x < b_threshold[2]) {
    pred = b_threshold[1] * x
  } else {
    pred = b_threshold[2] * b_threshold[1]
  }
  return(pred)
}

calc_mi_y <- function(po2, temp, b_threshold, Eo, n, Ao) {
  # calculate metabolic index
  ndata <- length(po2)
  kelvin = 273.15
  kb = 0.000086173324
  tref <- 15 + 273.15
  tempK <- temp + 273.15

  # create inv temp
  invtemp <- (1 / kb) * (1 / (tempK) - 1 / (tref))
  mi <- po2 * Ao * M ^ n * exp(Eo *invtemp)

  # apply threshold to scaled MI to get mu, and then sample from Tweedie to get observation
  y <- mu <- rep(NA, ndata)

  for (i in 1:ndata) {
    mu[i] <- exp(linear_threshold(mi[i], b_threshold))
    y[i] <- rtweedie(1,
                     mu = mu[i],
                     phi = phi,
                     power = s1)
  }
  return(list(mi = mi, y = y, mu = mu, invtemp = invtemp))
}


gen.data <- calc_mi_y(thedata$po2, thedata$temp, b_threshold, Eo, n, Ao)

thedata$mi <- gen.data$mi
thedata$y <- gen.data$y

# View simulated data and underlying  model
plot(thedata$mi, thedata$y)
points(thedata$mi, gen.data$mu, pch = 21, bg = "black")

### END CREATE DATA ###

### TMB FITTING  ####

### Set up prior covariances and means
priors <- function() {
  # load posterior means
  priormeans <- readRDS(file = "data/MI_posteriormeans.RDS")
  priorsigma <- readRDS(file = "data/MI_posteriorsigma.RDS")
  return(list(priormeans= priormeans, priorsigma = priorsigma))
}

prior <- priors()
thedata$tempK <- thedata$temp + 273.15
compile("src/baby_sdmTMB.cpp")
dyn.load(dynlib("src/baby_sdmTMB"))
data <- list(y = thedata$y,
             pO2 = thedata$po2,
             temp = thedata$tempK,
             tref = 15 + 273.15,
             W = 1000,
             ndata = ndata,
            priorsigma = prior$priorsigma,
             priormeans = prior$priormeans,
            nprior = 4L)

params <- list(Eo = Eo,
               n = n,
               logAomu = prior$priormeans[3],
               logAosigma = prior$priormeans[4],
               logAo = log(400),
               ln_phi = 3,
               thetaf = 0,
               b_threshold = c(-.5, 5)
)
obj <- MakeADFun(data = data, parameters = params, DLL = "baby_sdmTMB")
opt <- nlminb(obj$par, obj$fn, obj$gr)
rep <- sdreport(obj)
transformed <- summary(rep, "report")
print(transformed)
print(summary(rep, "fixed"))
sd_report <- TMB::sdreport(obj, getJointPrecision = get_joint_precision)

### Check Fits ####

Eofit <- summary(rep,"fixed")[1,1]
nfit <- summary(rep,"fixed")[2,1]
Aofit <- exp(summary(rep,"fixed")[5,1])
s2fit <- summary(rep, "fixed")[8,1]
scutfit <- summary(rep, "fixed")[9,1]

gen.data <- calc_mi_y(thedata$po2, thedata$temp, b_threshold= c(s2fit, scutfit), Eofit, nfit, Aofit)

points(thedata$mi, gen.data$mu, type = "p", pch = 21, bg = "red")
plot(thedata$mi, gen.data$mi)

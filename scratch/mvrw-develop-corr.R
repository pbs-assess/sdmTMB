library(ggplot2)
library(sdmTMB)

# as mvrw-develp but with unstructured correlation matrix rather than ordinal

set.seed(12928)
rho <- 0.6
stateDim <- 10
timeSteps <- 200
sds <- rlnorm(stateDim, meanlog = log(0.4), sdlog = 0.4)
sdObs <- rep(0.2, stateDim)
corrMat <- matrix(0.0, stateDim, stateDim)
offDiagLength <- (stateDim^2 - stateDim) / 2
rhoVec <- runif(offDiagLength, min = 0.1, max = 0.9)
corrMat[lower.tri(corrMat)] <- rhoVec
corrMat[lower.tri(corrMat)] <- t(corrMat)[lower.tri(corrMat)]
corrMat[diag(corrMat)] <- 1
# for (i in 1:stateDim) {
#   for (j in 1:stateDim) {
#     corrMat[i, j] <- rho^abs(i - j)
#   }
# }
SigmaRaw <- corrMat * (sds %o% sds)
#ensure Sigma pos definite
Sigma  <- Matrix::nearPD(SigmaRaw, conv.tol = 1e-7)$mat
d <- matrix(NA, timeSteps, stateDim)
obs <- d
d[1, ] <- rnorm(stateDim, 0, sds) # initial state
i <- 1
obs[i, ] <- d[i, ] + rnorm(stateDim, rep(0, stateDim), sdObs)
for (i in 2:timeSteps) {
  d[i, ] <- d[i - 1, ] + MASS::mvrnorm(1, rep(0, stateDim), Sigma = Sigma)
  obs[i, ] <- d[i, ] + rnorm(stateDim, rep(0, stateDim), sdObs)
}
matplot(d, type = "l")
truth <- d
matpoints(obs)

d <- data.frame(
  y = reshape2::melt(obs)[,3],
  year = rep(1:nrow(obs), stateDim),
  group = rep(letters[1:stateDim], each = timeSteps)
)
head(d)

dir.create("scratch/version-checks/", showWarnings = FALSE)

## Fits

rstudioapi::restartSession()
f <- "scratch/version-checks/fit-0.4.3.9002.rds"
if (!file.exists(f)) {
  remotes::install_github("sdmTMB/sdmTMB", ref = "96e5a92c") # before gengamma_Q
  library(sdmTMB)
  mesh <- make_mesh(pcod, c("X", "Y"), cutoff = 15)
  pcod$os <- rep(log(0.01), nrow(pcod)) # offset
  fit <- sdmTMB(
    data = pcod,
    formula = density ~ s(depth_scaled, k = 3),
    mesh = mesh,
    offset = pcod$os,
    family = tweedie(link = "log"),
    time = "year",
    time_varying = ~ 1,
    time_varying_type = 'ar1',
    extra_time = c(2006, 2008, 2010, 2012, 2014, 2016),
    spatiotemporal = "off"
  )
  saveRDS(fit, f)
}

rstudioapi::restartSession()
f <- "scratch/version-checks/fit-0.4.3.9003.rds"
if (!file.exists(f)) {
  remotes::install_github("sdmTMB/sdmTMB", ref = "031116ee")
  library(sdmTMB)
  mesh <- make_mesh(pcod, c("X", "Y"), cutoff = 15)
  pcod$os <- rep(log(0.01), nrow(pcod)) # offset
  fit <- sdmTMB(
    data = pcod,
    formula = density ~ s(depth_scaled, k = 3),
    mesh = mesh,
    offset = pcod$os,
    family = tweedie(link = "log"),
    time = "year",
    time_varying = ~ 1,
    time_varying_type = 'ar1',
    extra_time = c(2006, 2008, 2010, 2012, 2014, 2016),
    spatiotemporal = "off"
  )
  fit$version <- package_version("0.4.3.9003") # mistake!! entered as 0.4.2.9003 temporarily
  saveRDS(fit, f)
}

## Tests

devtools::install(".")
rstudioapi::restartSession()
library(sdmTMB)
library(testthat)

check_version <- function(f) {
  fit <- readRDS(f)
  print(fit$version)

  cli::cli_progress_message("Running update_version()")
  fit2 <- sdmTMB:::update_version(fit)
  pcod$os <- rep(log(0.01), nrow(pcod)) # offset

  cli::cli_progress_message("Running predict()")
  p <- predict(fit2, newdata = pcod)
  expect_equal(nrow(p), 2143)

  cli::cli_progress_message("Running residuals()")
  r <- residuals(fit2)
  expect_equal(length(r), 2143)

  cli::cli_progress_message("Running get_index()")
  p <- predict(fit2, newdata = subset(pcod, year == 2011), return_tmb_object = TRUE)
  ind <- get_index(p)
  expect_true(sum(is.na(ind$se)) == 0L)
}

f <- "scratch/version-checks/fit-0.4.3.9002.rds"
check_version(f)

f <- "scratch/version-checks/fit-0.4.3.9003.rds"
check_version(f)

## current
mesh <- make_mesh(pcod, c("X", "Y"), cutoff = 15)
pcod$os <- rep(log(0.01), nrow(pcod)) # offset
f <- "scratch/version-checks/fit-main.rds"
fit <- sdmTMB(
  data = pcod,
  formula = density ~ s(depth_scaled, k = 3),
  mesh = mesh,
  offset = pcod$os,
  family = tweedie(link = "log"),
  time = "year",
  time_varying = ~ 1,
  time_varying_type = 'ar1',
  extra_time = c(2006, 2008, 2010, 2012, 2014, 2016),
  spatiotemporal = "off"
)
saveRDS(fit, f)
check_version(f)

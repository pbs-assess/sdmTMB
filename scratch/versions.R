library(testthat)
dir.create("scratch/version-checks/", showWarnings = FALSE)

## Fits

f <- "scratch/version-checks/fit-0.4.3.9003.rds"
if (!file.exists(f)) {
  remotes::install_github("pbs-assess/sdmTMB", ref = "031116ee")
  mesh <- make_mesh(pcod, c("X", "Y"), cutoff = 15)
  pcod$os <- rep(log(0.01), nrow(pcod)) # offset
  fit <- sdmTMB(
    data = pcod,
    formula = density ~ s(depth_scaled, k = 3),
    mesh = mesh,
    offset = pcod$os,
    family = tweedie(link = "log"),
    time = "year",
    extra_time = c(2006, 2008, 2010, 2012, 2014, 2016),
    spatiotemporal = "ar1"
  )
  saveRDS(fit, f)
}

## Tests

devtools::load_all()
f <- "scratch/version-checks/fit-0.4.3.9003.rds"
fit <- readRDS(f)
fit2 <- update_version(fit)
pcod$os <- rep(log(0.01), nrow(pcod)) # offset
p <- predict(fit2, newdata = pcod)
expect_equal(nrow(p), 2143)

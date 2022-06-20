# comparing VAST and sdmTMB
#To do:
  #Distributions: Tweedie, delta-gamma, delta-lognormal
  #Spatiotemporal structures: IID, AR1 RW
    #VAST:2=RW; 3=constant among years as fixed effect; 4=AR1
  #To compare: Parameter estimates, predictions, COG, index values (w/ and w/out bias correction)

test_that("VAST logit-link delta-gamma", {
  ##VAST setup
  install.packages("devtools")
  library(devtools)
  library(VAST)
  ##sdmTMB setup
  install.packages("sdmTMB")
  library(sdmTMB)
  ##Add effort as variable (1 when using CPUE instead of observed weight as the response)
  pcod$effort <- 1
  if (FALSE) { # to check offset/effort
    set.seed(1)
    pcod$effort <- exp(rnorm(nrow(pcod), mean = 0, sd = 0.2))
  }
  ##VAST model
  #Make grid (=mesh)
  qcs_grid_ll <- qcs_grid
  qcs_grid_ll$Y <- qcs_grid_ll$Y * 1000
  qcs_grid_ll$X <- qcs_grid_ll$X * 1000
  qcs_grid_ll <- subset(qcs_grid_ll, year == min(qcs_grid_ll$year))
  input_grid <- cbind(Lat = qcs_grid_ll$Y, Lon = qcs_grid_ll$X, Area_km2 = 4)
  # with sp:
  coordinates(qcs_grid_ll) <- ~ X + Y
  proj4string(qcs_grid_ll) <- CRS("+proj=utm +zone=9")
  qcs_grid_ll <- as.data.frame(spTransform(qcs_grid_ll, CRS("+proj=longlat +datum=WGS84")))
  #make settings
  FieldConfig <- matrix(c("IID", "IID", "IID", "IID", "IID", "IID"),
                        ncol = 2, nrow = 3, dimnames = list(c("Omega", "Epsilon", "Beta"), c("Component_1", "Component_2")))
  RhoConfig <- c("Beta1" = 0, "Beta2" = 0, "Epsilon1" = 0, "Epsilon2" = 0)
  settings <- make_settings(
    n_x = 185, # number of vertices in the SPDE mesh
    Region = "User",
    purpose = "index2", # index of abundance with Gamma for positive catches
    fine_scale = TRUE, # use bilinear interpolation from the INLA 'A' matrix
    zone = 9,
    FieldConfig = FieldConfig,
    RhoConfig = RhoConfig,
    ObsModel = c(2, 0),
    bias.correct = FALSE,
    use_anisotropy = FALSE,
    max_cells = Inf,
    knot_method = "grid"
  )
  #fit model
  fit <- fit_model(
    settings = settings,
    Lat_i = pcod[, "lat"],
    Lon_i = pcod[, "lon"],
    t_i = pcod[, "year"],
    b_i = pcod[, "density"],
    a_i = pcod[, "effort"],
    input_grid = input_grid,
    working_dir = paste0(here("doc", "appendix-VAST"), "/")
  )
  ##sdmTMB model
  #make mesh (=grid)
  mesh <- make_mesh(pcod, xy_cols = c("X", "Y"),
                    mesh = fit$spatial_list$MeshList$isotropic_mesh)
  #fit model
  fit_sdmTMB <- sdmTMB(
    density ~ 0 + as.factor(year),
    data=pcod,
    mesh= mesh,
    offset=log(pcod$effort),
    family=delta_gamma(),
    time="year",
    silent=FALSE,
    control=sdmTMBcontrol(newton_loops=1)
  )

  ###Metrics to compare
  ##AIC
  #VAST
  #part of fit object
  #sdmTMB
  extractAIC.sdmTMB <- function(fit, scale, k = 2, ...) {
    L <- logLik(fit)
    edf <- attr(L, "df")
    return(c(edf, c(-2 * L + k * edf)))
  }
  extractAIC.sdmTMB(fit_sdmTMB)
  expect_equal(fit$parameter_estimates$AIC, fit2$AIC, tolerance=)
  ##Index
  #VAST
  plot <- plot(fit)
  index_vast <- plots$Index$table
  #sdmTMB
  p <- predict(fit2, newdata = qcs_grid, return_tmb_object = TRUE)
  index_sdm <- get_index(p, bias_correct = FALSE, area = rep(4, nrow(qcs_grid)))
  expect_equal(index_sdm, index_vast, tolerance=)
  ##Parameter estimates (which parameters?)
  ##Predictions
  ##Center of gravity
  })





##
test_that("VAST tweedie", {
  ##VAST setup
  install.packages("devtools")
  library(devtools)
  library(VAST)
  ##sdmTMB setup
  install.packages("sdmTMB")
  library(sdmTMB)
  ##Add effort as variable (1 when using CPUE instead of observed weight as the response)
  pcod$effort <- 1
  if (FALSE) { # to check offset/effort
    set.seed(1)
    pcod$effort <- exp(rnorm(nrow(pcod), mean = 0, sd = 0.2))
  }
  ##VAST model
  #Make grid (=mesh)
  qcs_grid_ll <- qcs_grid
  qcs_grid_ll$Y <- qcs_grid_ll$Y * 1000
  qcs_grid_ll$X <- qcs_grid_ll$X * 1000
  qcs_grid_ll <- subset(qcs_grid_ll, year == min(qcs_grid_ll$year))
  input_grid <- cbind(Lat = qcs_grid_ll$Y, Lon = qcs_grid_ll$X, Area_km2 = 4)
  # with sp:
  coordinates(qcs_grid_ll) <- ~ X + Y
  proj4string(qcs_grid_ll) <- CRS("+proj=utm +zone=9")
  qcs_grid_ll <- as.data.frame(spTransform(qcs_grid_ll, CRS("+proj=longlat +datum=WGS84")))
  #make settings
  ##??? Why is this done this way? Can I make it all IID to keep it simple?
  FieldConfig <- matrix(c("0", "0", "IID", "Identity", "IID", "IID", "IID", "Identity"),
                        ncol = 2, nrow = 4,
                        dimnames = list(
                          c("Omega", "Epsilon", "Beta", "Epsilon_year"),
                          c("Component_1", "Component_2")
                        )
  )
  RhoConfig <- c("Beta1" = 3, "Beta2" = 0, "Epsilon1" = 0, "Epsilon2" = 0)
  settings <- make_settings(
    n_x = 205, # number of vertices in the SPDE mesh
    Region = "User",
    purpose = "index2", # use recommended defaults for an index of abundance
    fine_scale = TRUE, # use bilinear interpolation from the INLA 'A' matrix
    zone = 9,
    FieldConfig = FieldConfig,
    RhoConfig = RhoConfig,
    ObsModel = c(10, 2), # use the Tweedie distribution as the observation model
    bias.correct = FALSE,
    use_anisotropy = FALSE,
    max_cells = Inf, # use all grid cells from the extrapolation grid
    knot_method = "grid" # or "samples"
  )
  #fit model
  fit <- fit_model(
    settings = settings,
    Lat_i = pcod[, "lat"],
    Lon_i = pcod[, "lon"],
    t_i = pcod[, "year"],
    b_i = pcod[, "density"],
    a_i = pcod[, "effort"],
    input_grid = input_grid,
  )
  ##sdmTMB model
  #make mesh (=grid)
  mesh <- make_mesh(pcod, xy_cols = c("X", "Y"),
                    mesh = fit$spatial_list$MeshList$isotropic_mesh)
  #fit model
  ##Change this to just
  fit_sdmTMB <- sdmTMB(density ~ 0 + as.factor(year),
                data = pcod, mesh = mesh, family = tweedie(link = "log"),
                time = "year", spatiotemporal = "iid", spatial = "on")

  ###Metrics to compare

})

###For AR1,
test_that("VAST tweedie AR1", {
  ##VAST setup
  install.packages("devtools")
  library(devtools)
  library(VAST)
  ##sdmTMB setup
  install.packages("sdmTMB")
  library(sdmTMB)
  ##Add effort as variable (1 when using CPUE instead of observed weight as the response)
  pcod$effort <- 1
  if (FALSE) { # to check offset/effort
    set.seed(1)
    pcod$effort <- exp(rnorm(nrow(pcod), mean = 0, sd = 0.2))
  }
  ##VAST model
  #Make grid (=mesh)
  qcs_grid_ll <- qcs_grid
  qcs_grid_ll$Y <- qcs_grid_ll$Y * 1000
  qcs_grid_ll$X <- qcs_grid_ll$X * 1000
  qcs_grid_ll <- subset(qcs_grid_ll, year == min(qcs_grid_ll$year))
  input_grid <- cbind(Lat = qcs_grid_ll$Y, Lon = qcs_grid_ll$X, Area_km2 = 4)
  # with sp:
  coordinates(qcs_grid_ll) <- ~ X + Y
  proj4string(qcs_grid_ll) <- CRS("+proj=utm +zone=9")
  qcs_grid_ll <- as.data.frame(spTransform(qcs_grid_ll, CRS("+proj=longlat +datum=WGS84")))
  #make settings
  FieldConfig <- matrix(c("0", "0", "IID", "Identity", "IID", "IID", "IID", "Identity"),
                        ncol = 2, nrow = 4,
                        dimnames = list(
                          c("Omega", "Epsilon", "Beta", "Epsilon_year"),
                          c("Component_1", "Component_2")
                        )
  )
  RhoConfig <- c("Beta1" = 0, "Beta2" = 0, "Epsilon1" = , "Epsilon2" = 0)
  settings <- make_settings(
    n_x = 205, # number of vertices in the SPDE mesh
    Region = "User",
    purpose = "index2", # use recommended defaults for an index of abundance
    fine_scale = TRUE, # use bilinear interpolation from the INLA 'A' matrix
    zone = 9,
    FieldConfig = FieldConfig,
    RhoConfig = RhoConfig,
    ObsModel = c(10, 2), # use the Tweedie distribution as the observation model
    bias.correct = FALSE,
    use_anisotropy = FALSE,
    max_cells = Inf, # use all grid cells from the extrapolation grid
    knot_method = "grid" # or "samples"
  )
  #fit model
  fit <- fit_model(
    settings = settings,
    Lat_i = pcod[, "lat"],
    Lon_i = pcod[, "lon"],
    t_i = pcod[, "year"],
    b_i = pcod[, "density"],
    a_i = pcod[, "effort"],
    input_grid = input_grid,
  )
  ##sdmTMB model
  #make mesh (=grid)
  mesh <- make_mesh(pcod, xy_cols = c("X", "Y"),
                    mesh = fit$spatial_list$MeshList$isotropic_mesh)
  #fit model
  fit_sdmTMB <- sdmTMB(density ~ 0 + as.factor(year),
                 data = pcod, mesh = mesh, family = tweedie(link = "log"),
                 time = "year", spatiotemporal = "iid", spatial = "on")

###metrics to compare
})

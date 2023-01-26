
# test_that("barrier mesh", {
#   skip_on_cran()
#   skip_on_ci()
#   skip_if_not_installed("INLA")
#   skip_if_not_installed("sf")
#   skip_if_not_installed("dplyr")
#
#
#   library(dplyr)
#   library(sf)
#   crs_utm9 <- 3156 # Pick a projection, here UTM9
#
#   suppressWarnings(
#   sf::st_crs(bc_coast) <- 4326 # 'WGS84'; necessary on some installs
#   )
#   bc_coast <- sf::st_transform(bc_coast, crs_utm9)
#
#   # Project our survey data coordinates:
#   survey <- pcod_2011 %>% select(lon, lat, density) %>%
#     st_as_sf(crs = 4326, coords = c("lon", "lat")) %>%
#     st_transform(crs_utm9)
#
#   # Plot our coast and survey data:
#   if (require("ggplot2", quietly = TRUE)) {
#   ggplot(bc_coast) +
#     geom_sf() +
#     geom_sf(data = survey, size = 0.5)
#   }
#   # Note that a barrier mesh won't don't much here for this
#   # example data set, but we nonetheless use it as an example.
#
#   surv_utm_coords <- st_coordinates(survey)
#
#   # Then we will scale coordinates to km so the range parameter
#   # is on a reasonable scale for estimation:
#   pcod_2011$X1000 <- surv_utm_coords[,1] / 1000
#   pcod_2011$Y1000 <- surv_utm_coords[,2] / 1000
#
#   spde <- make_mesh(pcod_2011, xy_cols = c("X1000", "Y1000"),
#                     cutoff = 10)
#   plot(spde)
#
#   # Add on the barrier mesh component:
#   bspde <- add_barrier_mesh(
#     spde, bc_coast, range_fraction = 0.1,
#     proj_scaling = 1000, plot = TRUE
#   )
#
#   # In the above, the grey dots are the centre of triangles that are in the
#   # ocean. The red crosses are centres of triangles that are over land. The
#   # spatial range will be assumed to be 0.1 (`range_fraction`) over land compared
#   # to over water.
#
#   # We can make a more advanced plot if we want:
#   mesh_df_water <- bspde$mesh_sf[bspde$normal_triangles, ]
#   mesh_df_land <- bspde$mesh_sf[bspde$barrier_triangles, ]
#
#   if (require("ggplot2", quietly = TRUE)) {
#   ggplot(bc_coast) +
#     geom_sf() +
#     geom_sf(data = mesh_df_water, size = 1, colour = "blue") +
#     geom_sf(data = mesh_df_land, size = 1, colour = "green")
#   }
#
#   fit <- sdmTMB(density ~ 1, data = pcod_2011, mesh = spde,
#                 family = tweedie(link = "log"))
#   fit
#
#   fit2 <- sdmTMB(density ~ 1, data = pcod_2011, mesh = bspde,
#                 family = tweedie(link = "log"))
#   fit2
#
#   (t1 <- tidy(fit, effects = "ran_pars", conf.int = TRUE))
#   (t2 <- tidy(fit2, effects = "ran_pars", conf.int = TRUE))
#
#   # if land is limiting range and in the regular spde, than range should increase with the barrier
#   expect_gt(t2$estimate[t2$term=="range"], t1$estimate[t1$term=="range"])
#   # not sure we'd always expect sigma_O to increase? But at least it has changed
#   expect_gt(t2$estimate[t2$term=="sigma_O"], t1$estimate[t1$term=="sigma_O"])
#
#   # try with time
#   fit <- sdmTMB(density ~ 1, data = pcod_2011, mesh = spde,
#                 time = "year", share_range = FALSE,
#                 family = tweedie(link = "log"))
#   fit
#
#   fit2 <- sdmTMB(density ~ 1, data = pcod_2011, mesh = bspde,
#                  time = "year", share_range = FALSE,
#                  family = tweedie(link = "log"))
#   fit2
#
#   (t1 <- tidy(fit, effects = "ran_pars", conf.int = TRUE))
#   (t2 <- tidy(fit2, effects = "ran_pars", conf.int = TRUE))
#
#   # if land is limiting range and in the regular spde, than range should increase with the barrier
#   expect_gt(t2$estimate[t2$term=="range"][1], t1$estimate[t1$term=="range"][1])
#   expect_gt(t2$estimate[t2$term=="range"][2], t1$estimate[t1$term=="range"][2])
#   # not sure we'd always expect sigma_O to increase? But at least it has changed
#   expect_gt(t2$estimate[t2$term=="sigma_O"], t1$estimate[t1$term=="sigma_O"])
#   # not sure what we would expect sigma_E to do?
#   # expect_??(t2$estimate[t2$term=="sigma_E"], t1$estimate[t1$term=="sigma_E"])
# })

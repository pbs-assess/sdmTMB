# d <- readRDS("~/src/gfsynopsis/report/data-cache/pacific-cod.rds")
# # d <- readRDS("~/src/gfsynopsis/report/data-cache/pacific-ocean-perch.rds")
# # d <- readRDS("~/src/gfsynopsis/report/data-cache/silvergray-rockfish.rds")
# # d <- readRDS("~/src/gfsynopsis/report/data-cache/redbanded-rockfish.rds")
# d <- d$survey_sets
# # d <- dplyr::filter(d, survey_series_id == 1)
# # d <- gfplot::get_survey_sets("pacific cod", ssid = 1)
#
#
# library(dplyr)
# # dat <- gfplot:::tidy_survey_sets(d, "SYN WCVI", years = seq(1, 1e6))
# dat <- gfplot:::tidy_survey_sets(d, "SYN QCS", years = seq(1, 1e6))
# # dat <- gfplot:::tidy_survey_sets(dat, "HBLL OUT S", years = seq(1, 3000), density_column = "density_ppkm2")
# # d <- d %>% select(year, longitude, latitude, depth_m, density_kgpm2)
#
# dat <- mutate(dat, density = density*1000*1000)
# dat <- filter(dat, !is.na(depth))
# dat <- gfplot:::scale_survey_predictors(dat)
# dat <- select(dat, -X10, -Y10)
#
# pcod <- dat
# usethis::use_data(pcod, internal = FALSE, overwrite = TRUE)
#
# grid_locs <- gfplot:::make_prediction_grid(filter(dat, year == 2003), survey = "SYN QCS", cell_width = 2)$grid
# # grid_locs <- gfplot:::make_prediction_grid(filter(dat, year == 2004), survey = "SYN WCVI", cell_width = 2)$grid
# grid_locs <- rename(grid_locs, depth = akima_depth)
# grid_locs$year <- NULL
# qcs_grid <- grid_locs
#
# # # Expand the prediction grid to create a slice for each time:
# # original_time <- sort(unique(dat$year))
# # nd <- do.call("rbind",
# #   replicate(length(original_time), qcs_grid, simplify = FALSE))
# # nd[["year"]] <- rep(original_time, each = nrow(qcs_grid))
#
# # qcs_grid <- nd
# usethis::use_data(qcs_grid, internal = FALSE, overwrite = TRUE)
#
# pcod_2011 <- subset(pcod, year >= 2011)
# pcod_mesh_2011 <- sdmTMB::make_mesh(pcod_2011, xy_cols = c("X", "Y"), cutoff = 20)
# usethis::use_data(pcod_mesh_2011, internal = FALSE, overwrite = TRUE)
#
# usethis::use_data(pcod_2011, internal = FALSE, overwrite = TRUE)
#
# # pcod_spde <- make_spde(pcod$X, pcod$Y, n_knots = 140)
# # plot_spde(pcod_spde)
# # m <- sdmTMB(
# #  pcod, density ~ 0 + as.factor(year) + depth_scaled + depth_scaled2,
# #  time = "year", spde = pcod_spde, family = tweedie(link = "log"), silent = FALSE,
# #   anisotropy = TRUE
# # )
# # plot_anisotropy(m)
# #
# # pcod$resids <- residuals(m) # randomized quantile residuals
# # hist(pcod$resids)
# # qqnorm(pcod$resids);abline(a = 0, b = 1)
# #
# # library(ggplot2)
# # ggplot(pcod, aes(X, Y, col = resids)) + scale_colour_gradient2() +
# #   geom_point() + facet_wrap(~year)
# #
# #
# # predictions <- predict(m, newdata = grid_locs)
# #
# # plot_map <- function(dat, column) {
# #   ggplot(dat, aes_string("X", "Y", fill = column)) +
# #     geom_raster() +
# #     facet_wrap(~year) +
# #     coord_fixed()
# # }
# #
# # plot_map(predictions$data, "exp(est)") +
# #   scale_fill_viridis_c(trans = "sqrt") +
# #   ggtitle("Prediction (fixed effects + all random effects)")
# #
# # plot_map(predictions$data, "exp(est_fe)") +
# #   ggtitle("Prediction (fixed effects only)") +
# #   scale_fill_viridis_c(trans = "sqrt")
# #
# # plot_map(predictions$data, "est_re_s") +
# #   ggtitle("Spatial random effects only") +
# #   scale_fill_gradient2()
# #
# # plot_map(predictions$data, "est_re_st") +
# #   ggtitle("Spatiotemporal random effects only") +
# #   scale_fill_gradient2()
# #
# # ind <- get_index(predictions, bias_correct = FALSE) # not bias correcting for speed
# #
# # library(ggplot2)
# # scale <- 2*2/(1000*1000)
# # ggplot(ind, aes(year, est*scale)) + geom_line() +
# #   geom_ribbon(aes(ymin = lwr*scale, ymax = upr*scale), alpha = 0.4)
# #
# # knitr::kable(ind, format = "pandoc")
# #
# #
#
# dd <- readRDS("../gfsynopsis-2021/report/data-cache-feb-2023/pacific-cod.rds")$survey_sets |>
#   dplyr::filter(survey_abbrev %in% "SYN QCS")
# dd$area_swept1 <- dd$doorspread_m * dd$tow_length_m
# dd$area_swept2 <- dd$doorspread_m * (dd$speed_mpm * dd$duration_min)
# dd$area_swept <- ifelse(!is.na(dd$area_swept1), dd$area_swept1, dd$area_swept2)
# new <- dplyr::select(dd, year, longitude, latitude, catch_weight, density_kgpm2, area_swept) |>
#   mutate(present = ifelse(catch_weight > 0, 1L, 0L))
#
# # old <- sdmTMB::pcod
#
# # new2 <- sdmTMB::add_utm_columns(new) |>
#   # filter(year <= 2017)
#   # rename(X2 = X, Y2 = Y)
#
# hist(new2$X)
# hist(old$X)
# range(old$year)
# range(new2$year)
#
# old$X_round <- round(old$X, digits = 2L)
# old$Y_round <- round(old$Y, digits = 2L)
#
# new2$X_round <- round(new2$X, digits = 2L)
# new2$Y_round <- round(new2$Y, digits = 2L)
#
# nrow(old)
# nrow(new)
#
# x <- left_join(old, select(new2, -X, -Y), by = c("X_round", "Y_round"))
# sum(is.na(x$area_swept))
#
#
# ##########

# library(dplyr)
# dd <- readRDS("../gfsynopsis-2021/report/data-cache-feb-2023/pacific-cod.rds")$survey_sets |>
#   dplyr::filter(survey_abbrev %in% "SYN QCS")
# saveRDS(dd, file = "data-raw/pcod-raw-no-git.rds")
# dd$area_swept1 <- dd$doorspread_m * dd$tow_length_m
# dd$area_swept2 <- dd$doorspread_m * (dd$speed_mpm * dd$duration_min)
# dd$area_swept <- ifelse(!is.na(dd$area_swept1), dd$area_swept1, dd$area_swept2)
# d <- dplyr::select(dd, year, longitude, latitude, catch_weight, density_kgpm2, area_swept, depth_m) |>
#   mutate(present = ifelse(catch_weight > 0, 1L, 0L)) |>
#   rename(depth = depth_m) |>
#   filter(!is.na(depth)) # 1 row
# d <- sdmTMB::add_utm_columns(d)
# # plot(d$catch_weight / d$area_swept, d$density_kgpm2)
# mean_log_depth <- mean(log(d$depth), na.rm = TRUE)
# sd_log_depth <- sd(log(d$depth), na.rm = TRUE)
# d <- mutate(d, depth_scaled = (log(depth) - mean_log_depth) / sd_log_depth)
# sum(is.na(d$depth))
# sum(is.na(d$depth_scaled))
# glimpse(d)
#
# d <- mutate(d, density = density_kgpm2 * 1e3 * 1e3)
# d <- mutate(d, area_swept = area_swept / 1e3 / 1e3)
# d <- select(d, year, longitude, latitude, X, Y, present, catch_weight, area_swept, density, depth, depth_scaled)
# pcod2 <- d
# glimpse(pcod2)
# digits <- 5L
# pcod2 <- pcod2 |> mutate(
#   area_swept = round(area_swept, digits),
#   density = round(density, digits),
#   longitude = round(longitude, digits),
#   latitude = round(latitude, digits),
#   X = round(X, digits),
#   Y = round(Y, digits),
#   depth_scaled = round(depth_scaled, digits)
# )
# pcod <- pcod2
# pcod_2011 <- dplyr::filter(pcod2, year >= 2011) # faster fitting
# usethis::use_data(pcod, internal = FALSE, overwrite = TRUE)
# usethis::use_data(pcod_2011, internal = FALSE, overwrite = TRUE)
#
# pcod_mesh_2011 <- sdmTMB::make_mesh(pcod_2011, xy_cols = c("X", "Y"), cutoff = 20)
# usethis::use_data(pcod_mesh_2011, internal = FALSE, overwrite = TRUE)
#
# qcs_grid <- gfplot::synoptic_grid |>
#   filter(survey == "SYN QCS") |>
#   dplyr::select(X, Y, depth)
#
# qcs_grid <- qcs_grid |> mutate(
#   X = round(X, digits),
#   Y = round(Y, digits),
#   depth = round(depth, digits)
# )
#
# usethis::use_data(qcs_grid, internal = FALSE, overwrite = TRUE)

## something non pcod??

library(dplyr)

make_pkg_data_trawl <- function(spp, survey) {
  dd <- readRDS(paste0("../gfsynopsis-2021/report/data-cache-feb-2023/", spp, ".rds"))$survey_sets |>
    dplyr::filter(survey_abbrev %in% survey)
  saveRDS(dd, file = paste0("data-raw/", spp, "-", gsub(" ", "-", survey), "-raw-no-git.rds"))
  dd$area_swept1 <- dd$doorspread_m * dd$tow_length_m
  dd$area_swept2 <- dd$doorspread_m * (dd$speed_mpm * dd$duration_min)
  dd$area_swept <- ifelse(!is.na(dd$area_swept1), dd$area_swept1, dd$area_swept2)
  d <- dplyr::select(dd, year, longitude, latitude, catch_weight, density_kgpm2, area_swept, depth_m) |>
    mutate(present = ifelse(catch_weight > 0, 1L, 0L)) |>
    rename(depth = depth_m) |>
    filter(!is.na(depth)) # 1 row
  d <- sdmTMB::add_utm_columns(d)
  d <- mutate(d, area_swept = area_swept / 1e3 / 1e3)
  d <- select(d, year, longitude, latitude, X, Y, present, catch_weight, area_swept, depth)
  digits <- 5L
  d <- d |> mutate(
    area_swept = round(area_swept, digits),
    longitude = round(longitude, digits),
    latitude = round(latitude, digits),
    X = round(X, digits),
    Y = round(Y, digits)
  )
  as.data.frame(d)
}
dogfish <- make_pkg_data_trawl("north-pacific-spiny-dogfish", "SYN WCVI")
usethis::use_data(dogfish, internal = FALSE, overwrite = TRUE)

pcod_mesh_2011 <- sdmTMB::make_mesh(pcod_2011, xy_cols = c("X", "Y"), cutoff = 20)
usethis::use_data(pcod_mesh_2011, internal = FALSE, overwrite = TRUE)

# digits <- 5L
# qcs_grid <- gfplot::synoptic_grid |>
#   filter(survey == "SYN QCS") |>
#   dplyr::select(X, Y, depth)
# qcs_grid <- qcs_grid |>
#   mutate(
#     X = round(X, digits),
#     Y = round(Y, digits),
#     depth = round(depth, 0)
#   ) |>
#   distinct()
# usethis::use_data(qcs_grid, internal = FALSE, overwrite = TRUE)

wcvi_grid <- gfplot::synoptic_grid |>
  filter(survey == "SYN WCVI") |>
  dplyr::select(X, Y, depth)
wcvi_grid <- wcvi_grid |>
  mutate(
    X = round(X, digits),
    Y = round(Y, digits),
    depth = round(depth, 0)
  ) |>
  distinct()
usethis::use_data(wcvi_grid, internal = FALSE, overwrite = TRUE)

glimpse(wcvi_grid)

library(sdmTMB)
mesh <- make_mesh(dogfish, c("X", "Y"), cutoff = 10)
fit <- sdmTMB(
  catch_weight ~ s(depth),
  data = dogfish,
  anisotropy = TRUE,
  mesh = mesh,
  time = "year",
  share_range = FALSE,
  offset = log(dogfish$area_swept),
  silent = FALSE,
  family = tweedie()
)
sanity(fit)
summary(fit)
plot_anisotropy(fit)

# -------------- longline?

spp <- "yelloweye-rockfish"
survey <- "HBLL OUT S"
dd <- readRDS(paste0("../gfsynopsis-2021/report/data-cache-feb-2023/", spp, ".rds"))$survey_sets |>
  dplyr::filter(survey_abbrev %in% survey)
glimpse(dd)
bait <- readRDS("../gfsynopsis-2021/report/data-cache-july-2023/bait-counts.rds")

nrow(dd)
dd <- left_join(dd, bait)
nrow(dd)
sum(is.na(dd$count_bait_only))
sum(dd$count_bait_only > dd$hook_count)
saveRDS(dd, file = paste0("data-raw/", spp, "-", gsub(" ", "-", survey), "-raw-no-git.rds"))

d <- dplyr::select(dd, year, longitude, latitude, catch_count, hook_count, depth_m, bait_returned = count_bait_only) |>
  rename(depth = depth_m) |>
  filter(!is.na(depth))
d <- sdmTMB::add_utm_columns(d)
d <- select(d, year, longitude, latitude, X, Y, catch_count, hook_count, bait_returned, depth)
digits <- 5L
d <- d |> mutate(
  catch_count = as.integer(catch_count),
  hook_count = as.integer(hook_count),
  bait_returned = as.integer(bait_returned),
  longitude = round(longitude, digits),
  latitude = round(latitude, digits),
  X = round(X, digits),
  Y = round(Y, digits),
  depth = round(depth, 1)
)
yelloweye <- as.data.frame(d)
usethis::use_data(yelloweye, internal = FALSE, overwrite = TRUE)
glimpse(yelloweye)

library(sdmTMB)
mesh <- make_mesh(yelloweye, c("X", "Y"), cutoff = 10)
fit <- sdmTMB(
  catch_count ~ s(depth),
  data = yelloweye,
  anisotropy = TRUE,
  share_range = FALSE,
  mesh = mesh,
  time = "year",
  offset = log(yelloweye$hook_count),
  silent = FALSE,
  family = nbinom2()
)
sanity(fit)
summary(fit)
plot_anisotropy(fit)

hbll_s_grid <- gfplot::hbll_s_grid$grid
hbll_s_grid <- hbll_s_grid |>
  mutate(
    longitude = X,
    latitude = Y,
    depth = round(depth, 1)
  ) |>
  select(-X, -Y) |>
  distinct()

hbll_s_grid <- sdmTMB::add_utm_columns(hbll_s_grid)
hbll_s_grid <- hbll_s_grid |> select(X, Y, depth)
hbll_s_grid <- hbll_s_grid |>
  mutate(
    X = round(X, 5),
    Y = round(Y, 5)
  )
usethis::use_data(hbll_s_grid, internal = FALSE, overwrite = TRUE)

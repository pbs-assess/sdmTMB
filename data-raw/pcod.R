d <- readRDS("~/src/gfsynopsis/report/data-cache/pacific-cod.rds")
d <- d$survey_sets
# d <- dplyr::filter(d, survey_series_id == 1)
# d <- gfplot::get_survey_sets("pacific cod", ssid = 1)


library(dplyr)
dat <- gfplot:::tidy_survey_sets(d, "SYN QCS", years = seq(1, 1e6))
# dat <- gfplot:::tidy_survey_sets(dat, "HBLL OUT S", years = seq(1, 3000), density_column = "density_ppkm2")
# d <- d %>% select(year, longitude, latitude, depth_m, density_kgpm2)

dat <- mutate(dat, density = density*1000*1000)
dat <- filter(dat, !is.na(depth))
dat <- gfplot:::scale_survey_predictors(dat)
dat <- select(dat, -X10, -Y10)

pcod <- dat
usethis::use_data(pcod, internal = FALSE, overwrite = TRUE)

grid_locs <- gfplot:::make_prediction_grid(filter(dat, year == 2003), survey = "SYN QCS", cell_width = 2)$grid
grid_locs <- rename(grid_locs, depth = akima_depth)
grid_locs$year <- NULL
qcs_grid <- grid_locs
usethis::use_data(qcs_grid, internal = FALSE, overwrite = TRUE)

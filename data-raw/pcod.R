d <- readRDS("~/src/gfsynopsis/report/data-cache/pacific-cod.rds")
d <- d$survey_sets
d <- dplyr::filter(d, survey_series_id == 1)
# d <- gfplot::get_survey_sets("pacific cod", ssid = 1)

library(dplyr)
d <- d %>% select(year, longitude, latitude, depth_m, density_kgpm2)

pcod <- d
usethis::use_data(pcod, internal = FALSE, overwrite = TRUE)

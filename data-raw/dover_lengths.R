library(dplyr)
d <- readRDS("../gfsynopsis-2023/report/data-cache-2024-05/dover-sole.rds")$survey_samples
d_sets <- readRDS("../gfsynopsis-2023/report/data-cache-2024-05/dover-sole.rds")$survey_sets
d <- filter(d, survey_abbrev == "SYN WCVI", !is.na(length))
d <- left_join(distinct(d), distinct(select(d_sets, survey_abbrev, year, longitude, latitude, fishing_event_id)))
d <- select(d, year, length, lon = longitude, lat = latitude, trawl_id = fishing_event_id)
d <- filter(d, !is.na(lon), !is.na(lat))
# saveRDS(d, "scratch/dover_lengths.rds")

dover_lengths <- d
usethis::use_data(dover_lengths, overwrite = TRUE)

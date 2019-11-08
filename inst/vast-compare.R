library(sp)
library(ggplot2)
library(dplyr)
library(sdmTMB)

grid <- readr::read_csv("~/Downloads/EBSThorsonGrid.csv")
coordinates(grid) <- c("Lon", "Lat")
proj4string(grid) <- CRS("+proj=longlat +datum=WGS84")
crs <- paste0(
  "+proj=aea +lat_1=55 +lat_2=65 +lat_0=50 +lon_0=-154 +x_0=0",
  " +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs"
)
grid <- as.data.frame(spTransform(grid, CRS(crs)))
colnames(grid) <- tolower(colnames(grid))
grid$lat <- grid$lat / 10000
grid$lon <- grid$lon / 10000
plot(grid$lon, grid$lat, pch = ".")

d <- readRDS("~/Downloads/EBSpollock.rds")
coordinates(d) <- c("lon", "lat")
proj4string(d) <- CRS("+proj=longlat +datum=WGS84")
d <- as.data.frame(spTransform(d, CRS(crs)))
d$lat <- d$lat / 10000
d$lon <- d$lon / 10000

spde <- make_spde(d$lon, d$lat, n_knots = 250)
plot_spde(spde)
m <- sdmTMB(
  data = d, formula = cpue_kg_km2 ~ 0 + as.factor(year), spde = spde,
  family = tweedie(link = "log"), time = "year", silent = FALSE
)
saveRDS(m, file = "~/Downloads/pollock-ebs-sdmTMB.rds")
# m$tmb_obj$retape()

original_time <- sort(unique(d$year))
nd <- do.call(
  "rbind",
  replicate(length(original_time), grid, simplify = FALSE)
)
nd[["year"]] <- rep(original_time, each = nrow(grid))
nd <- dplyr::select(nd, year, lon, lat)

p <- predict(m, newdata = nd, return_tmb_object = TRUE, xy_cols = c("lon", "lat"))

g <- ggplot(p$data, aes(lon, lat)) + geom_tile(aes(fill = est), width = 0.375, height = 0.375) +
  facet_wrap(vars(year)) +
  scale_fill_viridis_c() +
  coord_fixed()
ggsave("~/Desktop/sdmTMB-Tweedie-pollock.png", width = 14, height = 14)

ind <- get_index(p, bias_correct = FALSE)
saveRDS(ind, file = "~/Downloads/sdmtmb-pollock-index.rds")

comparison <- readr::read_csv("~/Downloads/21740_estimate_summary.csv")

comp <- comparison %>%
  select(year = Year, est = VAST_mt, se = VAST_CV) %>%
  mutate(lwr = exp(log(est) - 1.96 * se)) %>%
  mutate(upr = exp(log(est) + 1.96 * se)) %>%
  mutate(type = "VAST: Poisson link")

geo_mean_ratio <- exp(mean(log(ind$est))) / exp(mean(log(comp$est)))
.ind <- mutate(ind, type = "sdmTMB: Tweedie") %>%
  mutate(
    est = est / geo_mean_ratio, lwr = lwr / geo_mean_ratio,
    upr = upr / geo_mean_ratio
  )
comp <- .ind %>% bind_rows(comp)

ggplot(comp, aes(year, y = est, ymin = lwr, ymax = upr)) +
  geom_line(aes(color = type)) +
  geom_ribbon(aes(fill = type), alpha = 0.2) +
  theme_light() +
  xlab("") + ylab("Estimate") + labs(colour = "Type", fill = "Type")
ggsave("~/Desktop/VAST-Poisson-link-vs-sdmTMB-Tweedie.pdf", width = 8, height = 5)

fit_sdmTMB_westcoast <- function(species_rds, survey,
  species_name = "", n_knots = 200, cell_width = 2,
  anisotropy = FALSE, silent = TRUE, bias_correct = FALSE) {

  d <- readRDS(species_rds)$survey_sets
  d <- dplyr::filter(d, !(year == 2014 & survey_abbrev == "SYN WCHG")) # not used
  col <- if (grepl("SYN", survey)) "density_kgpm2" else "density_ppkm2"
  dat <- gfplot:::tidy_survey_sets(d, survey, years = seq(1, 1e6),
    density_column = col)

  if (mean(dat$present) < 0.05) error("Not enough data.")

  .scale <- if (grepl("SYN", survey)) 1000 else 1 # for computational stability
  dat <- dplyr::mutate(dat, density = density * .scale)
  if (any(is.na(dat$depth)))
    dat <- gfplot:::interp_survey_bathymetry(dat)$data
  # dat <- dplyr::filter(dat, !is.na(depth))
  dat <- gfplot:::scale_survey_predictors(dat)

  if (grepl("SYN", survey)) {
    grid_locs <- gfplot:::make_prediction_grid(
      dplyr::filter(dat, year == max(dat$year)), survey = survey,
      cell_width = cell_width)$grid
    grid_locs <- dplyr::rename(grid_locs, depth = akima_depth)
  } else {
    grid_locs <- if (surv == "HBLL OUT N") gfplot::hbll_n_grid$grid else gfplot::hbll_s_grid$grid
  }
  grid_locs$year <- NULL

  spde <- make_spde(dat$X, dat$Y, n_knots = n_knots)
  m <- sdmTMB(
    formula = density ~ 0 + as.factor(year), #+ depth_scaled + depth_scaled2,
    data = dat, time = "year", spde = spde, family = tweedie(link = "log"),
    anisotropy = anisotropy, silent = silent)
  predictions <- predict(m, newdata = grid_locs)

  index <- get_index(predictions, bias_correct = bias_correct)

  # scale the biomass by the area and adjust units;
  # this will be cleaned up and integrated
  # scale <- 2 * 2 / 1000 # 2 x 2 km grid and converted from kg to tonnes
  index <- dplyr::mutate(index,
    cv = sqrt(exp(se^2) - 1)
  )

  list(
    data = dat,
    model = m,
    spde = spde,
    predictions = predictions,
    index = index,
    scale = .scale,
    survey = survey,
    species_name = species_name
  )
}

load_all("../gfsynopsis/")
spp <- gfsynopsis::get_spp_names()
spp <- dplyr::pull(dplyr::filter(spp, type %in% c("A", "B")), spp_w_hyphens)
# spp <- "yelloweye-rockfish"
survs <- c('SYN QCS', 'SYN HS', 'SYN WCHG', 'SYN WCVI')
# survs <- c('HBLL OUT N', 'HBLL OUT S')
# survs <- c('OTHER HS MSA')
#survs <- c('SYN QCS')
# survs <- c('HS MSA')
all <- expand.grid(spp = spp, survs = survs,
  stringsAsFactors = FALSE)
library(sdmTMB)
library(foreach)
cores <- min(nrow(all), parallel::detectCores())
cl <- parallel::makeCluster(cores)
doParallel::registerDoParallel(cl)
out <- foreach::foreach(sp = all$spp, surv = all$survs,
  .packages = c("gfplot", "sdmTMB")) %dopar% {
    tryCatch(fit_sdmTMB_westcoast(
      here::here(paste0(
        "../gfsynopsis/report/data-cache/", sp, ".rds")),
      species_name = sp,
      survey = surv, n_knots = 200, bias_correct = FALSE,
      anisotropy = FALSE
    ), error = function(e) NA)
  }
doParallel::stopImplicitCluster()
saveRDS(out, file = "inst/spt-index-out-2018-10-19.rds")

library(ggplot2)
library(dplyr)

index <- purrr::map_df(out, function(x) {
  if (length(x) > 1L)
    data.frame(x$index, species = x$species_name, survey = x$survey,
      stringsAsFactors = FALSE) %>% tibble::as.tibble()
})
saveRDS(index, file = "inst/spt-index-out-no-depth.rds")

index$survey <- factor(index$survey,
  levels = c('SYN WCHG', 'SYN HS', 'SYN QCS', 'SYN WCVI', 'HBLL OUT N', 'HBLL OUT S'))
ggplot(index, aes(year, est)) + geom_line() +
  geom_ribbon(aes(ymin = lwr, ymax = upr), alpha = 0.4) +
  xlab('Year') + ylab('Biomass estimate (metric tonnes)') +
  facet_grid(species~survey, scales = "free")

design_based <- purrr::map_df(unique(index$species), function(sp) {
  message(sp)
  .d <- readRDS(here::here(paste0(
    "../gfsynopsis/report/data-cache/", sp, ".rds")))
  .d$survey_index
})

index <- index %>%
  group_by(survey, species) %>%
  mutate(
    lwr = lwr / exp(mean(log(est))),
    upr = upr / exp(mean(log(est))),
    est = est / exp(mean(log(est)))
  ) %>%
  ungroup()

des <- design_based %>%
  group_by(survey_abbrev, species_common_name) %>%
  mutate(
    lowerci = lowerci / exp(mean(log(biomass))),
    upperci = upperci / exp(mean(log(biomass))),
    biomass = biomass / exp(mean(log(biomass)))
  ) %>%
  # mutate(
  #   lowerci = lowerci / (1000*1000),
  #   upperci = upperci / (1000*1000),
  #   biomass = biomass / (1000*1000)
  # ) %>%
  ungroup() %>%
  select(year, biomass, lowerci, upperci, survey_abbrev, species_common_name, re) %>%
  filter(survey_abbrev %in% unique(index$survey)) %>%
  rename(est = biomass, lwr = lowerci, upr = upperci, survey = survey_abbrev, species = species_common_name, cv = re) %>%
  mutate(species = gsub(" ", "-", species)) %>%
  mutate(species = gsub("/", "-", species)) %>%
  mutate(type = "Design based")

index$type <- "Spatiotemporal"
ind <- suppressWarnings(bind_rows(index, des))
inds <- group_by(ind, survey, species, type) %>%
  summarise(
    max_cv = max(cv, na.rm = TRUE) < 1,
    max_est = max(est) < 50,
    cv_na = all(!is.na(cv)))
inds <- inds %>% filter(max_cv, max_est, cv_na)
ind <- semi_join(ind, inds)
ind <- filter(ind, species != "pacific-hake")

# ind <- filter(ind, cv < 1)
# ind <- filter(ind, est < 20)
# # ind <- filter(ind, species != "big-skate")
# ind <- filter(ind, !(year == 2014 & survey == "SYN WCHG"))
ind$survey <- factor(ind$survey,
  levels = c('SYN WCHG', 'SYN HS', 'SYN QCS', 'SYN WCVI', 'HBLL OUT N', 'HBLL OUT S'))
saveRDS(ind, file = "inst/index-estimates-2018-10-19.rds")

g <- ggplot(ind, aes(year, est, fill = type)) +
  geom_line(aes(colour = type)) +
  geom_point(aes(colour = type)) +
  geom_ribbon(aes(ymin = lwr, ymax = upr), alpha = 0.4) +
  xlab('Year') + ylab('Relative biomass estimate') +
  facet_grid(species~survey, scales = "free_y") +
  scale_x_continuous(breaks = seq(2000, 2020, 5)) +
  labs(colour = "Type", fill = "Type") #+
  # gfplot::theme_pbs()

ggsave("inst/surv-2018-10-19-no-depth-200-knots.pdf", width = 9.5, height = 65, limitsize = FALSE)


# plot_spde(out$spde)
#
# plot_anisotropy(out$model)
#
# out$data$resids <- residuals(out$model) # randomized quantile residuals
# hist(out$data$resids)
# qqnorm(out$data$resids);abline(a = 0, b = 1)
# ggplot(out$data, aes(X, Y, col = resids)) + scale_colour_gradient2() +
#   geom_point() + facet_wrap(~year) + coord_fixed()
#
# plot_map <- function(dat, column) {
#   ggplot(dat, aes_string("X", "Y", fill = column)) +
#     geom_raster() +
#     facet_wrap(~year) +
#     coord_fixed()
# }
#
# plot_map(out$predictions$data, "exp(est)") +
#   scale_fill_viridis_c(trans = "sqrt") +
#   ggtitle("Prediction (fixed effects + all random effects)")
#
# plot_map(out$predictions$data, "exp(est_fe)") +
#   ggtitle("Prediction (fixed effects only)") +
#   scale_fill_viridis_c(trans = "sqrt")
#
# plot_map(out$predictions$data, "est_re_s") +
#   ggtitle("Spatial random effects only") +
#   scale_fill_gradient2()
#
# plot_map(out$predictions$data, "est_re_st") +
#   ggtitle("Spatiotemporal random effects only") +
#   scale_fill_gradient2()
#
#

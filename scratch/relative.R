library(sdmTMB) # 'scale' branch
library(ggplot2)
library(dplyr)

if (Sys.info()[["user"]] == "seananderson") {
  # dat <- readRDS("../gfsynopsis-2021/report/data-cache-april-2022/longnose-skate.rds")$survey_sets
  dat <- readRDS("../gfsynopsis-2021/report/data-cache-april-2022/yellowtail-rockfish.rds")$survey_sets
  dat <- filter(dat, survey_abbrev == "SYN QCS")
  dat$density <- dat$density_kgpm2
  dat <- add_utm_columns(dat)
} else {
  dat <- pcod
}

mesh <- make_mesh(dat, c("X", "Y"), cutoff = 12)
# plot(mesh)

fit <- sdmTMB(
  density ~ 0 + as.factor(year),
  data = dat,
  time = "year",
  mesh = mesh,
  family = tweedie(link = "log"),
  spatial = "on",
  spatiotemporal = "off", # fast testing
  silent = FALSE,
  control = sdmTMBcontrol(newton_loops = 1L)
)

nd <- replicate_df(qcs_grid, "year", unique(dat$year))
p <- predict(fit, newdata = nd, return_tmb_object = TRUE)

ind_abs <- get_index(p, area = 4, bias_correct = FALSE)
ind_rel_mean <- get_index(p, type = "relative-mean", area = 4, bias_correct = FALSE)
ind_rel_one <- get_index(p, type = "relative-year1", area = 4, bias_correct = FALSE)
ind <- rbind(
  mutate(ind_abs, type = "absolute"),
  mutate(ind_rel_mean, type = "relative-mean"),
  mutate(ind_rel_one, type = "relative-year1")
)
ind_labs <- group_by(ind, type) |>
  summarise(
    lab_x = min(year) + 0.5,
    lab_y = max(upr),
    mean_se = paste0("Mean SE = ", round(mean(se), 2))
  )

ggplot(ind, aes(year, est, fill = type)) +
  geom_pointrange(aes(ymin = lwr, ymax = upr), position = position_dodge(width = 0.5)) +
  facet_wrap(vars(type), scales = "free_y") +
  geom_text(aes(lab_x, lab_y, label = mean_se), data = ind_labs, colour = "black", hjust = 0) +
  theme_light()

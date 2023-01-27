library(sdmTMB) # 'scale' branch
library(ggplot2)
library(dplyr)

if (Sys.info()[['user']] == "seananderson") {
  # dat <- readRDS("../gfsynopsis-2021/report/data-cache-april-2022/longnose-skate.rds")$survey_sets
  dat <- readRDS("../gfsynopsis-2021/report/data-cache-april-2022/yellowtail-rockfish.rds")$survey_sets
  dat <- filter(dat, survey_abbrev == "SYN QCS")
  dat$density <- dat$density_kgpm2
  dat <- add_utm_columns(dat)
} else {
  dat <- pcod
}

mesh <- make_mesh(dat, c("X", "Y"), cutoff = 12)
plot(mesh)

fit <- sdmTMB(
  density ~ 0 + as.factor(year),
  data = dat,
  time = "year",
  mesh = mesh,
  family = tweedie(link = "log"),
  spatial = "on",
  spatiotemporal = "off",
  silent = FALSE
)

nd <- replicate_df(qcs_grid, "year", unique(dat$year))
p <- predict(fit, newdata = nd, return_tmb_object = TRUE)

system
ind_abs <- get_index(p, relative = FALSE, area = 4, bias_correct = FALSE)
ind_rel <- get_index(p, relative = TRUE, area = 4, bias_correct = FALSE)
ind <- rbind(
  mutate(ind_abs, type = "absolute"),
  mutate(ind_rel, type = "relative")
)

ggplot(ind, aes(year, est, fill = type, colour = type)) +
  geom_line() +
  geom_ribbon(aes(ymin = lwr, ymax = upr), alpha = 0.4, colour = NA) +
  facet_wrap(vars(type), scales = "free_y")

group_by(ind, type) |>
  summarise(mean_se = mean(se))

library(sdmTMB) # 'scale' branch
library(ggplot2)
library(dplyr)

mesh <- make_mesh(pcod, c("X", "Y"), cutoff = 15)
plot(mesh)

fit <- sdmTMB(
  density ~ 0 + as.factor(year),
  data = pcod,
  time = "year",
  mesh = mesh,
  family = tweedie(link = "log"),
  spatial = "on",
  spatiotemporal = "iid",
  silent = FALSE
)

nd <- replicate_df(qcs_grid, "year", unique(pcod$year))
p <- predict(fit, newdata = nd, return_tmb_object = TRUE)

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

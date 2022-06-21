library(sdmTMB)
mesh <- make_mesh(pcod, c("X", "Y"), cutoff = 20)

# typical modular approach:
fit <- sdmTMB(
  density ~ s(depth),
  data = pcod, mesh = mesh,
  time = "year",
  family = tweedie(link = "log")
)
pred <- predict(fit, newdata = qcs_grid, return_tmb_object = TRUE)
index <- get_index(pred, area = 4)

# do index calculations during fitting like in VAST:
fit2 <- sdmTMB(
  density ~ s(depth),
  data = pcod, mesh = mesh,
  time = "year",
  family = tweedie(link = "log"),
  do_index = TRUE, #<
  index_args = list(area = 4), #< or a vector of length nrow(qcs_grid)
  predict_args = list(newdata = qcs_grid) #<
)
index2 <- get_index(fit2)

index
index2

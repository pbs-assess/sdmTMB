library(sdmTMB)
library(ggplot2)

mesh <- make_mesh(pcod_2011, c("X", "Y"), cutoff = 10)
options(sdmTMB.cores = 4)
tictoc::tic()
fit <- sdmTMB(density ~ 0 + depth_scaled + depth_scaled2 + as.factor(year),
  data = pcod_2011, time = "year", mesh = mesh, family = tweedie(link = "log"))
tictoc::toc()

# ------

tests <- expand.grid(iter = 1, i = 1:6)
df <- purrr::pmap_dfr(tests, function(iter, i) {
  cat(i, "\n")
  options(sdmTMB.cores = i)
  tictoc::tic()
  fit <- sdmTMB(density ~ 0 + depth_scaled + depth_scaled2 + as.factor(year),
    data = pcod_2011, time = "year", mesh = mesh, family = tweedie(link = "log"))
  a <- tictoc::toc()
  data.frame(iter = iter, n_threads = i, time = unname(a$toc - a$tic))
})
ggplot(df, aes(n_threads, time)) + geom_point()

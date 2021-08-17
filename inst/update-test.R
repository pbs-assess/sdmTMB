old <- readRDS("~/Downloads/mod-imm-biomass-pacific-cod-tv-depth-only-1n3n4n16.rds")
load_all()
# new <- readRDS("~/Downloads/mod-imm-biomass-pacific-cod-tv-depth-only-1n3n4n16-2021.rds")
new <- readRDS("~/Downloads/mod-imm-biomass-pacific-cod-tv-depth-only-1n3n4n16-2021b.rds")

old2 <- sdmTMB::update_model(old, xy_cols = c("X", "Y"))

class(new) <- "sdmTMB"
print(new)
print(old2)

# # for (i in seq_along(names(new$tmb_data))) {
#   i <- 0
#   i <- i + 1
#   this <- names(new$tmb_data)[i]
#   print(this)
#   if (!identical(new$tmb_data[[this]], old2$tmb_data[[this]])) {
#     print(new$tmb_data[[this]])
#     print(old2$tmb_data[[this]])
#     print(this)
#     if (i > length(names(new$tmb_data))) stop("STOP!")
#   }
# # }

nrow(new$data)
nrow(old$data)

pold <- predict(old2)
class(new) <- "sdmTMB"
pnew <- predict(new)
plot(pold$est, pnew$est)

pold <- predict(old2, newdata = NULL)
class(new) <- "sdmTMB"
pnew <- predict(new, newdata = NULL)
plot(pold$est, pnew$est)

new$tmb_params$ln_tau_E

d <- subset(pcod, year >= 2015)
pcod_spde <- make_mesh(d, c("X", "Y"), cutoff = 100)
new <- sdmTMB(density ~ 0 + as.factor(year),
  data = d, time = "year", spde = pcod_spde, family = tweedie(link = "log"),
  do_fit = T)

od <- names(old$tmb_data)
nd <- names(new$tmb_data)

missing <- nd[!nd %in% od]
missing

od <- names(old$tmb_params)
nd <- names(new$tmb_params)

missing <- nd[!nd %in% od]
missing

od <- names(old$tmb_map)
nd <- names(new$tmb_map)

missing <- nd[!nd %in% od]
missing

od <- names(old$tmb_random)
nd <- names(new$tmb_random)

missing <- nd[!nd %in% od]
missing

old2 <- sdmTMB::update_model(old, xy_cols = c("X", "Y"))

p <- predict(old2, newdata = NULL)
p <- predict(old2)
tidy(old2)
print(old2)



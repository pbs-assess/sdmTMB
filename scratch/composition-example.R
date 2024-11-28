library(dplyr)
library(ggplot2)

d <- readr::read_csv("~/Downloads/Sablefish_StRS_Data/StRS_LWMSO.csv")

d[1:4,] |> glimpse()

nrow(d)
d <- filter(d, !is.na(SPECIMEN_AGE))
nrow(d)

glimpse(d)
unique(d$SPECIES_CODE)

names(d) <- tolower(names(d))

glimpse(d)
unique(d$gear_code)

d <- select(d, year, lon = slon, lat = slat, age = specimen_age, catch_weight, fishing_event_id)

glimpse(d)
d$year <- as.integer(d$year)

ggplot(d, aes(lon, lat)) + geom_point(aes(colour = age)) +
  facet_wrap(~year) +
  scale_colour_viridis_c()

table(d$year, d$age)

cuts <- seq(0, 120, 2)

d$age_binned <- cuts[findInterval(d$age, cuts)]
d$age_binned[d$age_binned >= 40] <- 40

ggplot(d, aes(lon, lat)) + geom_point(aes(colour = age)) +
  facet_wrap(~year) +
  scale_colour_viridis_c()

ggplot(d, aes(lon, lat)) + geom_point(aes(colour = age)) +
  facet_grid(age_binned~year) +
  scale_colour_viridis_c()

dg <- group_by(d, year, fishing_event_id, lon, lat, age_binned) |>
  summarise(n = n(), .groups = "drop")

dg <- sdmTMB::add_utm_columns(dg, c("lon", "lat"))

ggplot(dg, aes(X, Y)) + geom_point(aes(colour = n)) +
  facet_grid(age_binned~year) +
  scale_colour_viridis_c()

table(dg$year, dg$age_binned)

dtest <- subset(dg, age_binned == 12)

library(sdmTMB)

mesh <- make_mesh(dtest, c("X", "Y"), cutoff = 10)
plot(mesh)
mesh$mesh$n

dtest$id <- factor(1:nrow(dtest))

fit <- sdmTMB(n ~ 0 + factor(year), mesh = mesh, data = dtest, family = poisson(),
  time = "year", silent = FALSE, spatiotemporal = "off",
  priors = sdmTMBpriors(matern_s = pc_matern(range_gt = 2, sigma_lt = 1)))
summary(fit)
sanity(fit)

sort(unique(dg$age_binned))

dtest <- subset(dg, age_binned %in% seq(4, 40, 2))
mesh <- make_mesh(dtest, c("X", "Y"), cutoff = 10)
dtest$age_binned_factor <- factor(dtest$age_binned)
nages <- length(unique(dtest$age_binned_factor))

table(dtest$year, dtest$age_binned_factor)

fit <- sdmTMB(
  n ~ 0 + factor(year) * age_binned_factor,
  mesh = mesh,
  data = dtest,
  family = poisson(),
  time = "year",
  silent = FALSE,
  spatial = "off",
  spatial_varying = ~ 0 + age_binned_factor,
  # control = sdmTMBcontrol(map = list(ln_tau_Z = factor(rep(1L, nages))), profile = TRUE),
  control = sdmTMBcontrol(map = list(ln_tau_Z = factor(c(1, 2, rep(3, nages - 2)))), profile = TRUE),
  spatiotemporal = "off",
  priors = sdmTMBpriors(matern_s = pc_matern(range_gt = 5, sigma_lt = 1))
)
summary(fit)
sanity(fit)

# get abundance density indexes for each year-age combination:

# pretend WCVI for now:
nd <- replicate_df(sdmTMB::wcvi_grid, "year", unique(dtest$year))
ages <- unique(dtest$age_binned_factor)
ages

# loop over ages:
ind_list <- lapply(ages, \(a) {
  print(a)
  nd$age_binned_factor <- a
  pred <- predict(fit, newdata = nd, return_tmb_object = TRUE)
  ind <- get_index(pred, area = 4)
  data.frame(ind, age_bin = a)
})
ind <- do.call(rbind, ind_list)

ggplot(ind, aes(year, est)) +
  geom_ribbon(aes(ymin = lwr, ymax = upr)) +
  geom_line(colour = "white") +
  facet_wrap(~age_bin, scales = "fixed")

ind$I_ct <- ind$est
ind <- group_by(ind, year) |> mutate(I_t = sum(I_ct))
ind <- mutate(ind, P_ct = I_ct / I_t)
ind <- group_by(ind, age_bin) |>
  mutate(sum_var = sum(se_natural^2))

nyrs <- length(unique(dtest$year))
nages <- length(unique(dtest$age_binned_factor))
x <- t(matrix(ind$est, nrow = nyrs, ncol = nages))
Index_ctl <- array(x, dim = c(nages, nyrs, 1))
x <- t(matrix(ind$se_natural, nrow = nyrs, ncol = nages))
SE_Index_ctl <- array(x, dim = c(nages, nyrs, 1))
# Calculate proportions, and total biomass
Prop_ctl = Index_ctl / outer(rep(1,dim(Index_ctl)[1]), apply(Index_ctl,MARGIN=2:3,FUN=sum))
Index_tl = apply(Index_ctl,MARGIN=2:3,FUN=sum)
SE_Index_tl = sqrt(apply(SE_Index_ctl^2,MARGIN=2:3,FUN=sum,na.rm=TRUE))
# Approximate variance for proportions, and effective sample size
Neff_ctl = var_Prop_ctl = array(NA,dim=dim(Prop_ctl))
for( cI in 1:dim(var_Prop_ctl)[1]){
  for( tI in 1:dim(var_Prop_ctl)[2]){
    for( lI in 1:dim(var_Prop_ctl)[3]){
      var_Prop_ctl[cI,tI,lI] = Index_ctl[cI,tI,lI]^2/Index_tl[tI,lI]^2 * (SE_Index_ctl[cI,tI,lI]^2/Index_ctl[cI,tI,lI]^2 - 2*SE_Index_ctl[cI,tI,lI]^2/(Index_ctl[cI,tI,lI]*Index_tl[tI,lI]) + SE_Index_tl[tI,lI]^2/Index_tl[tI,lI]^2 )
      var_Prop_ctl[cI,tI,lI] = ifelse( Index_ctl[cI,tI,lI]==0, 0, var_Prop_ctl[cI,tI,lI] )  # If dividing by zero, replace with 0
      # Covert to effective sample size
      Neff_ctl[cI,tI,lI] = Prop_ctl[cI,tI,lI] * (1-Prop_ctl[cI,tI,lI]) / var_Prop_ctl[cI,tI,lI]
    }}}
# Median effective sample size across categories
Neff_tl = apply(Neff_ctl, MARGIN=2:3, FUN=median, na.rm=TRUE)

# # var_Prop_ctl[cI,tI,lI] =
# #   Index_ctl[cI,tI,lI]^2/Index_tl[tI,lI]^2 *
# #   (SE_Index_ctl[cI,tI,lI]^2/Index_ctl[cI,tI,lI]^2 -
# #       2*SE_Index_ctl[cI,tI,lI]^2/(Index_ctl[cI,tI,lI]*Index_tl[tI,lI]) +
# #       SE_Index_tl[tI,lI]^2/Index_tl[tI,lI]^2 )
#
# ind <- ind |> mutate(
#   var_Pct = # Eq. 12
#     I_ct^2 / I_t^2 *
#     # Index_ctl[cI,tI,lI]^2/Index_tl[tI,lI]^2 *
#     (se_natural^2 / I_ct^2 -
#         # (SE_Index_ctl[cI,tI,lI]^2/Index_ctl[cI,tI,lI]^2 -
#         2 * (se_natural^2 / I_t * I_ct) +
#         # 2*SE_Index_ctl[cI,tI,lI]^2/(Index_ctl[cI,tI,lI]*Index_tl[tI,lI]) +
#         sum_var / I_t^2)
# )

# ind <- ind |> mutate(se_Pct = sqrt(var_Pct))

out <- reshape::melt(var_Prop_ctl) |> select(-X3) |> rename(age_bin_i = X1, year_i = X2) |>
  mutate(prop_se = sqrt(value)) |> select(-value)
out2 <- reshape::melt(Prop_ctl) |> select(-X3) |> rename(age_bin_i = X1, year_i = X2) |>
  mutate(prop = value) |> select(-value)
out3 <- reshape::melt(Neff_ctl) |> select(-X3) |> rename(age_bin_i = X1, year_i = X2) |>
  mutate(Neff = value) |> select(-value)
out4 <- reshape::melt(Neff_tl) |> select(-X2) |> rename(year_i = X1, Neff_median = value)
out <- left_join(out, out2) |> left_join(out3) |> left_join(out4)
# out

luy <- select(dtest, year) |> distinct() |> mutate(year_i = 1:n())
lua <- select(dtest, age_binned) |> distinct() |> mutate(age_bin_i = 1:n())

out <- left_join(out, luy) |> left_join(lua)

ggplot(out, aes(age_binned, prop)) +
  facet_wrap(~year) +
  geom_linerange(aes(ymin = prop - .5 * prop_se, ymax = prop + .5 * prop_se)) +
  geom_line() +
  coord_cartesian(ylim = c(0, NA))

group_by(out, year_i) |> summarise(neff_median = Neff_median[1])

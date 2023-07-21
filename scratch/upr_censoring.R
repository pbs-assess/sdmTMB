# Used in conjunction with a browser() here: https://github.com/joenomiddlename/Censored_Longline_RCode/blob/main/Case%20Study/Case_Study_Analysis_Standalone_After_Revisions.Rmd#L204
# Used this to get the test dataframe.
# i <- "yelloweye"
# sf::sf_use_s2(FALSE)
#
# censored_index_fun_sdmTMB(
#   data = Reduced_data_sp, survey_boundaries = survey_boundaries, species = i, M = 1000, return = T, ICR_adjust = F, cprop = 0.95, keep = T, use_upper_bound = T, upper_bound_quantile = 0.95, plot = F, allyears = F, station_effects = T, seed = 0, verbose = F, n_trajectories = 0,
#   preserve_inter_regional_differences = F, time_effect = "unstructured",
#   twentyhook_adjust = F
# )
#
# data$scale_fac <- scale_fac
# data$upper_bound <- upper_bound
#
# test_df <-
#   data %>%
#   dplyr::select(year, N_dat, prop_removed, obsHooksPerSet, scale_fac, upper_bound, low, high) %>%
#   dplyr::slice(1:40)

# write_csv(test_df, '../sdmTMB/scratch/upr_censoring_test_df.csv')

test_df <- read.csv(here::here("scratch", "upr_censoring_test_df.csv"))
# p_tk = proportion of baits removed in fishing event k of year t
# p_istar = true breakdown point for species i, which is cprop?


# # Equation S.4
# comp_factor_fun <-function(prop_hook, n_hooks) {
#     prop <- 1 - prop_hook
#     # if all hooks saturated - map to 1 hook
#     prop[which(prop == 0)] <- 1 / n_hooks[which(prop == 0)]
#     return(-log(prop) / (1 - prop))
#     }
#
# # using L176
# # Why are there cases with cprop = 1.1?
# get_scale_factor <- function(prop_removed, n_hooks, cprop) {
#   prop_hook = signif(((prop_removed - cprop) / (1 - cprop)), 5)
#   n_hooks = round((1 - cprop) * n_hooks)
#
#   scale_fac <- comp_factor_fun(prop_hook = prop_hook, n_hooks = n_hooks)
#
#   return(scale_fac)
# }

#' Calculate scaling factor for hook competition using censored method
#'
#' @param prop_removed A vector (numeric) containing the proportion of baits
#'   removed in each fishing event.
#' @param n_hooks A vector (integers) containing the number of hooks
#'   deployed in each fishing event.
#' @param pstar A single value (numeric) specifying the breakdown point of
#'   observed catch counts as a result of hook competition.
#' @return A vector containing the scale factor used to calculate the
#'   upper bound on catch counts of target species to improve convergence of
#'   censored method. See \doi{10.1139/cjfas-2022-0159}, including supplementary
#'   materials section S.3 for more details.
#' @noRd
get_scale_factor <- function(prop_removed, n_hooks, pstar) {
  # Apply pstar correction
  prop_hook <- signif(((prop_removed - pstar) / (1 - pstar)), 5)
  n_hooks <- round((1 - pstar) * n_hooks)
  # Calculate competition adjustment factor
  prop <- 1 - prop_hook
  prop[prop == 0] <- 1 / n_hooks[prop == 0] # if all hooks saturated - map to 1 hook
  -log(prop) / (1 - prop)
}

#' Calculate an upper bound on catch counts for the censored Poisson family
#'
#' @param prop_removed The proportion of baits removed in each fishing event
#'   from *any* species. I.e., the proportion of hooks returning without bait
#'   for any reason.
#' @param n_catch The observed catch counts on each fishing event of the target
#'   species.
#' @param n_hooks The number of hooks deployed on each fishing event
#' @param pstar A single value between `0 <= pstar <= 1` specifying the
#'   breakdown point of observed catch counts as a result of hook competition.
#'
#' @details `pstar` could be obtained via inspecting a GAM or other smoother fit
#'   with catch counts as the response an offset for log(hook count) and
#'   proportion of baits removed for each fishing event and checking when the
#'   curve drops off.
#'
#' @return A numeric vector of upper bound catch counts of the target species to
#'   improve convergence of censored method.
#'
#' @references See \doi{10.1139/cjfas-2022-0159} for more details.
#' @export
#'
#' @examples
dat <- structure(
  list(
    n_catch = c(78L, 63L, 15L, 6L, 7L, 11L, 37L, 99L, 34L, 100L, 77L, 79L,
      98L, 30L, 49L, 33L, 6L, 28L, 99L, 33L),
    prop_removed = c(
      0.61, 0.81, 0.96, 0.69, 0.99, 0.98, 0.25, 0.95, 0.89, 1, 0.95, 0.95,
      0.94, 1, 0.95, 1, 0.84, 0.3, 1, 0.99
    ), n_hooks = c(
      140L, 140L, 140L, 140L, 140L, 140L, 140L, 140L, 140L, 140L, 140L, 140L,
      140L, 140L, 140L, 140L, 140L, 140L, 140L, 140L
    )),
  class = "data.frame", row.names = c(NA, -20L)
)
upr <- get_upper_bound(dat$prop_removed, dat$n_catch, dat$n_hooks)
plot(dat$n_catch, upr, type = "n")
symbols(dat$n_catch, upr, circles = dat$prop_removed)
abline(0, 1)
plot(dat$prop_removed, upr)

get_upper_bound <- function(
    prop_removed,
    n_catch,
    n_hooks,
    pstar = 0.95) {

  assertthat::assert_that(
    is.numeric(prop_removed),
    is.numeric(n_catch),
    is.numeric(n_hooks)
  )
  assertthat::assert_that(length(prop_removed) == length(n_catch))
  assertthat::assert_that(length(prop_removed) == length(n_hooks))
  assertthat::assert_that(pstar >= 0)
  assertthat::assert_that(pstar <= 1)
  assertthat::assert_that(sum(is.na(c(prop_removed, n_catch, n_hooks))) == 0)

  removed_ind <- prop_removed > pstar
  # FIXME: add cli::cli_abort()?
  # probably need to throw error if length(removed_ind) == 0L
  upper_bound <- rep(0, length(prop_removed))
  scale_fac <- rep(0, length(prop_removed))

  # Calculate part of scaling factor used to get upper bound (part of eq. S.20)
  scale_fac[removed_ind] <-
    get_scale_factor(
      pstar = pstar,
      prop_removed = prop_removed[removed_ind],
      n_hooks = n_hooks[removed_ind]
    )

  # Upper bound for a species (part of eq. S.22)
  upper_bound[removed_ind] <- (prop_removed[removed_ind] - pstar) *
    n_hooks[removed_ind] * scale_fac[removed_ind]

  high <- n_catch
  high[prop_removed >= pstar] <- high[prop_removed >= pstar] +
    upper_bound[prop_removed >= pstar]
  round(high)
}

#
# # Testing:
# test_upper <-
#   get_upper_bound(
#     test_df$prop_removed,
#     test_df$N_dat,
#     test_df$obsHooksPerSet,
#     pstar = 0.95
#   )
#
# test_df$high
# test_upper$upr
#
# all.equal(test_df$high, test_upper$upr)

# cbind(test_df, test_upper$upr) %>% View()

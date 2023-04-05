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
#' @param prop_removed A vector (as numeric) containing the proportion of baits
#'   removed in fishing event
#' @param n_hooks An vector (as integer?) containing the number of hooks
#'   deployed in a fishing event
#' @param cprop A value (as numeric) specifying the breakdown point of observed
#'   catch counts as a result of hook competition.
#' @return An object (vector) containing the scale factor used to calculate the
#'   upper bound on catch counts of target species to improve convergence of
#'   censored method. See \doi{10.1139/cjfas-2022-0159}, including supplementary
#'   materials section S.3 for more details.
#' @noRd
get_scale_factor <- function(prop_removed, n_hooks, cprop) {
  # Apply cprop correction
  prop_hook <- signif(((prop_removed - cprop) / (1 - cprop)), 5)
  n_hooks <- round((1 - cprop) * n_hooks)

  # Calculate competition adjustment factor
  prop <- 1 - prop_hook
  prop[prop == 0] <- 1 / n_hooks[prop == 0] # if all hooks saturated - map to 1 hook
  -log(prop) / (1 - prop)
}


#' Calculate upper bound on catch counts for censored method
#'
#' @param prop_removed A numeric vector containing the proportion of baits
#'   removed in a fishing event.
#' @param n_catch A numeric vector containing observed catch counts.
#' @param n_hooks An numeric vector containing the number of hooks deployed in a
#'   fishing event
#' @param cprop A numeric value specifying the breakdown point of observed catch
#'   counts as a result of hook competition.
#'
#' @return A data frame with columns containing the `lwr` and `upr` bounds on
#'   catch counts of target species to improve convergence of censored method.
#'   See \doi{10.1139/cjfas-2022-0159} for more details.
#' @export
#'
#' @examples
#' dat <- structure(
#'   list(n_catch = c(
#'     0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L,
#'     0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 8L, 0L, 0L
#'   ), prop_removed = c(
#'     0.61, 0.81, 0.96, 0.69, -0.99, 0.98, -0.25, 0.95, 0.89, 1, 0.95, 0.95,
#'     0.94, 1, 0.95, 1, 0.84, 0.3, 1, 0.99
#'   ), n_hooks = c(
#'     140L, 140L, 140L, 140L, 140L, 140L, 140L, 140L, 140L, 140L, 140L, 140L,
#'     140L, 140L, 140L, 140L, 140L, 140L, 140L, 140L
#'   )),
#'   class = "data.frame", row.names = c(NA, -20L)
#' )
#' get_upper_bound(dat$prop_removed, dat$n_catch, dat$n_hooks)

get_upper_bound <- function(
    prop_removed,
    n_catch,
    n_hooks,
    cprop = 0.95) {

  assertthat::assert_that(is.numeric(prop_removed), is.numeric(n_catch), is.numeric(n_hooks))
  assertthat::assert_that(length(prop_removed) == length(n_catch))
  assertthat::assert_that(length(prop_removed) == length(n_hooks))
  assertthat::assert_that(cprop > 0)
  assertthat::assert_that(cprop < 1)
  assertthat::assert_that(sum(is.na(c(prop_removed, n_catch, n_hooks))) == 0)

  N_dat <- n_catch

  removed_ind <- prop_removed > cprop
  # FIXME: add cli::cli_abort()?
  # probably need to throw error if length(removed_ind) == 0L
  upper_bound <- rep(0, length(prop_removed))
  scale_fac <- rep(0, length(prop_removed))

  # Calculate part of scaling factor used to get upper bound (part of eq. S.20)
  scale_fac[removed_ind] <-
    get_scale_factor(cprop, prop_removed = prop_removed[removed_ind], n_hooks = n_hooks[removed_ind])

  # Upper bound for a species (part of eq. S.22)
  upper_bound[removed_ind] <- (prop_removed[removed_ind] - cprop) *
    n_hooks[removed_ind] * scale_fac[removed_ind]

  low <- N_dat
  high <- N_dat

  # FIXME: is there a reason that this is prop_removed > cprop, not prop_removed >= cprop?
  high[prop_removed >= cprop] <- high[prop_removed >= cprop] +
    upper_bound[prop_removed >= cprop]

  upper_bound <- round(upper_bound, 3)
  high <- round(high)
  data.frame(low = low, upr = high)
}

# Testing:
test_upper <-
  get_upper_bound(
    test_df$prop_removed,
    test_df$N_dat,
    test_df$obsHooksPerSet,
    cprop = 0.95
  )

test_df$high
test_upper$upr

all.equal(test_df$high, test_upper$upr)

# cbind(test_df, test_upper$upr) %>% View()

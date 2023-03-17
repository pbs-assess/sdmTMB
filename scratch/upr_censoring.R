# Used in conjunction with a browser() here: https://github.com/joenomiddlename/Censored_Longline_RCode/blob/main/Case%20Study/Case_Study_Analysis_Standalone_After_Revisions.Rmd#L204
# Used this to get the test dataframe.
i = "yelloweye"
sf::sf_use_s2(FALSE)

censored_index_fun_sdmTMB(
  data=Reduced_data_sp, survey_boundaries=survey_boundaries, species=i, M=1000, return=T, ICR_adjust=F, cprop=0.95, keep=T, use_upper_bound=T, upper_bound_quantile=0.95, plot=F, allyears=F, station_effects=T, seed=0, verbose=F, n_trajectories=0,
  preserve_inter_regional_differences = F, time_effect = 'unstructured',
  twentyhook_adjust = F
)

data$scale_fac <- scale_fac
data$upper_bound <- upper_bound

test_df <-
  data %>% dplyr::select(year, N_dat, prop_removed, obsHooksPerSet, scale_fac, upper_bound, low, high) %>%
  dplyr::slice(1:40)

#write_csv(test_df, '../sdmTMB/scratch/upr_censoring_test_df.csv')


test_df <- read.csv(here::here('scratch', 'upr_censoring_test_df.csv'))
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

#
#' Calculate scaling factor for hook competition using censored method
#'
#' @param prop_removed A vector (as numeric) containing the proportion of baits removed in fishing event
#' @param n_hooks An vector (as integer?) containing the number of hooks deployed in a fishing event
#' @param cprop A value (as numeric) specifying the breakdown point of observed catch counts as a result of hook competition.
#'
#' @return An object (vector) containing the scale factor used to calculate the upper bound on catch counts of target species to improve convergence of censored method. See \doi{10.1139/cjfas-2022-0159}, including supplementary materials section S.3 for more details.
#' @export # Should this just be an internal function?
#'
#' @examples
get_scale_factor <- function(prop_removed, n_hooks, cprop) {
  # Apply cprop correction
  prop_hook = signif(((prop_removed - cprop) / (1 - cprop)), 5)
  n_hooks = round((1 - cprop) * n_hooks)

  # Calculate competition adjustment factor
  prop <- 1 - prop_hook
  prop[which(prop == 0)] <- 1 / n_hooks[which(prop == 0)]  # if all hooks saturated - map to 1 hook
  scale_fac <- -log(prop) / (1 - prop)

  return(scale_fac)
}


#' Calculate upper bound on catch counts for censored method
#'
#' @param data A data frame.
#' @param prop_removed A column name (as character) of the column containing the proportion of baits removed in a fishing event.
#' @param N_dat A column name (as character) of the column containing observed catch counts.
#' @param n_hooks An column name (as character) containing the number of hooks deployed in a fishing event
#' @param cprop A value (as numeric) specifying the breakdown point of observed catch counts as a result of hook competition.
#'
#' @return A data frame with columns containing the `lwr` and `upr` bounds on catch counts of target species to improve convergence of censored method. See \doi{10.1139/cjfas-2022-0159} for more details.
#' @export
#'
#' @examples
get_upper_bound <- function(data, prop_removed = 'prop_removed', N_dat = 'N_dat', n_hooks = 'obsHooksPerSet', cprop = cprop) {
  prop_removed <- data[[prop_removed]]
  N_dat <- data[[N_dat]]
  n_hooks <- data[[n_hooks]]
  removed_ind <- prop_removed > cprop  # probably need to throw error if length(removed_ind) == 0L
  upper_bound <- rep(0, length(prop_removed))
  scale_fac <- rep(0, length(prop_removed))

  # Calculate part of scaling factor used to get upper bound (part of eq. S.20)
  scale_fac[removed_ind] <-
    get_scale_factor(cprop, prop_removed = prop_removed[removed_ind], n_hooks = n_hooks[removed_ind])

  # Upper bound for a species (part of eq. S.22)
  upper_bound[removed_ind] <- (prop_removed[removed_ind] - cprop) * n_hooks[removed_ind] * scale_fac[removed_ind]

  low  <- N_dat
  high <- N_dat

  high[which(prop_removed >= cprop)] <- high[which(prop_removed >= cprop)] + upper_bound[which(prop_removed >= cprop)]  # Is there a reason that this is prop_removed > cprop, not prop_removed >= cprop?

  scale_fac = round(scale_fac, 3)
  upper_bound <- round(upper_bound, 2)
  high <- round(high)

  censoring_bounds <- data.frame(low = low, upr = high)
  #test <- cbind(prop_removed, N_dat, n_hooks, cprop, removed_ind, upper_bound, scale_fac, low, high)
  return(censoring_bounds)
}

# Testing:
test_upper <-
  get_upper_bound(data = test_df, prop_removed = 'prop_removed', N_dat = 'N_dat',
                  n_hooks = 'obsHooksPerSet', cprop = 0.95)

test_df$high
test_upper$upr

all.equal(test_df$high, test_upper$upr)

cbind(test_df, test_upper$upr) %>% view

#' Process distribution column for observation-level families
#'
#' Internal helper function to handle distribution_column parameter in sdmTMB().
#' Creates family and link arrays for TMB when different observations have
#' different distributions. Supports mixed delta and non-delta families.
#' Uses unified approach for both single family and mixed distribution cases.
#'
#' @param family Either a single family object or named list of family objects
#' @param distribution_column Character string of column name in data, or NULL
#' @param data Data frame containing the distribution column
#' @param y_i Response matrix (n_obs x n_m)
#' @param n_m Number of model components (1 for single, 2 for delta)
#'
#' @return List containing:
#'   \item{tmb_family}{Family codes for TMB (array)}
#'   \item{tmb_link}{Link codes for TMB (array)}
#'   \item{tmb_poisson_link_delta}{Poisson link delta indicator (vector)}
#'   \item{d_i}{0-indexed vector mapping observations to phi parameter groups}
#'   \item{n_m_effective}{Effective number of model components}
#'
#' @keywords internal
#' @noRd
process_distribution_column <- function(family, distribution_column, data, y_i, n_m) {
  n_obs <- length(y_i)
  
  # Normalize inputs - convert single family to distribution_column format
  if (is.null(distribution_column)) {
    # Create fake distribution column and family list for single family
    family_list <- list("single" = family)
    dist_values <- rep("single", n_obs)
  } else {
    # Use provided distribution column
    if (!is.list(family) || is.null(names(family))) {
      cli::cli_abort("When `distribution_column` is provided, `family` must be a named list.")
    }
    family_list <- family
    dist_values <- data[[distribution_column]]
    
    # Validate distributions
    if (!all(dist_values %in% names(family_list))) {
      missing_dists <- unique(dist_values[!dist_values %in% names(family_list)])
      cli::cli_abort(paste("Distribution(s) not found in family names:", 
                          paste(missing_dists, collapse = ", ")))
    }
  }
  
  # Unified logic for all cases from here
  d_i <- as.integer(factor(dist_values)) - 1L  # 0-indexed for TMB
  obs_families <- family_list[dist_values]
  
  # Determine component requirements per observation
  is_delta <- vapply(obs_families, function(x) isTRUE(x$delta), logical(1))
  components_per_obs <- ifelse(is_delta, 2L, 1L)
  max_components <- max(components_per_obs)
  
  has_delta <- any(is_delta)
  has_non_delta <- any(!is_delta)
  
  if (has_delta && has_non_delta) {
    # Mixed case: use maximum components needed with proper usage tracking
    n_m_effective <- max_components  # Always 2 for mixed delta/non-delta
    if (!is.null(distribution_column)) {
      cli::cli_inform("Mixed delta and non-delta families detected. Using component-wise likelihood evaluation.")
    }
  } else if (has_delta) {
    # All delta: use n_m = 2
    n_m_effective <- 2L
  } else {
    # All non-delta: use provided n_m
    n_m_effective <- n_m
  }
  
  # Build arrays using vectorized logic
  family_array <- array(NA_integer_, dim = c(n_obs, n_m_effective))
  link_array <- array(NA_integer_, dim = c(n_obs, n_m_effective))
  component_usage <- array(FALSE, dim = c(n_obs, n_m_effective))
  
  # Create logical vector for delta vs non-delta observations
  is_delta_obs <- is_delta
  
  # Separate delta and non-delta indices for efficient processing
  delta_indices <- which(is_delta_obs)
  non_delta_indices <- which(!is_delta_obs)
  
  # Process non-delta observations (vectorized)
  if (length(non_delta_indices) > 0) {
    non_delta_families <- obs_families[non_delta_indices]
    family_array[non_delta_indices, 1] <- .valid_family[vapply(non_delta_families, function(x) x$family, character(1))]
    link_array[non_delta_indices, 1] <- .valid_link[vapply(non_delta_families, function(x) x$link, character(1))]
    component_usage[non_delta_indices, 1] <- TRUE
  }
  
  # Process delta observations (vectorized)
  if (length(delta_indices) > 0) {
    delta_families <- obs_families[delta_indices]
    # Extract first component (binomial part)
    family_array[delta_indices, 1] <- .valid_family[vapply(delta_families, function(x) x$family[1], character(1))]
    link_array[delta_indices, 1] <- .valid_link[vapply(delta_families, function(x) x$link[1], character(1))]
    component_usage[delta_indices, 1] <- TRUE
    
    # Extract second component (positive part) if needed
    if (n_m_effective > 1L) {
      family_array[delta_indices, 2] <- .valid_family[vapply(delta_families, function(x) x$family[2], character(1))]
      link_array[delta_indices, 2] <- .valid_link[vapply(delta_families, function(x) x$link[2], character(1))]
      component_usage[delta_indices, 2] <- TRUE
    }
  }
  
  # Vectorized poisson_link_delta extraction
  poisson_link_delta_vec <- as.integer(vapply(obs_families, 
    function(x) isTRUE(x$type == "poisson_link_delta"), logical(1)))
  
  return(list(
    tmb_family = family_array,
    tmb_link = link_array,
    tmb_poisson_link_delta = poisson_link_delta_vec,
    d_i = d_i,
    n_m_effective = n_m_effective,
    component_usage = component_usage,
    components_per_obs = components_per_obs
  ))
}

#' Safe family checking for both single and list family cases
#'
#' @param family Either a single family object or named list of family objects
#' @param distribution_column Character string of column name in data, or NULL
#' @return List with safe family information
#' @keywords internal
#' @noRd
safe_family_check <- function(family, distribution_column = NULL) {
  if (is.null(distribution_column)) {
    # Single family case
    return(list(
      first_family = family$family[1],
      all_families = family$family,
      is_single = TRUE,
      family_obj = family
    ))
  } else {
    # Multiple families case - cannot do simple checks
    return(list(
      first_family = NULL,
      all_families = NULL,
      is_single = FALSE,
      family_obj = family
    ))
  }
}

#' Validate early family constraints before distribution processing
#'
#' @param family Either a single family object or named list of family objects
#' @param distribution_column Character string of column name in data, or NULL
#' @param y_i Response matrix
#' @param delta Logical indicating delta model
#' @param upr Upper censoring bound for censored Poisson
#' @param experimental Experimental parameters list
#' @param control Control parameters
#' @param data Data frame (needed for row count validation)
#' @keywords internal
#' @noRd
validate_early_family_constraints <- function(family, distribution_column, y_i, delta, upr, experimental, control, data) {
  family_info <- safe_family_check(family, distribution_column)
  
  if (family_info$is_single) {
    # Only validate for single family case - mixed families handled later
    first_fam <- family_info$first_family
    
    # Censored Poisson validation
    if (first_fam == "censored_poisson") {
      if ("lwr" %in% names(experimental) || "upr" %in% names(experimental)) {
        cli_abort("Detected `lwr` or `upr` in `experimental`. `lwr` is no longer needed and `upr` is now specified as `control = sdmTMBcontrol(censored_upper = ...)`.")
      }
      if (is.null(control$censored_upper)) {
        cli_abort("`censored_upper` must be defined in `control = sdmTMBcontrol()` to use the censored Poisson distribution.")
      }
      assert_that(length(control$censored_upper) == nrow(data))
      assert_that(mean(control$censored_upper - y_i, na.rm = TRUE) >= 0)
    }
    
    # Gamma/lognormal positive value validation
    if (first_fam %in% c("Gamma", "lognormal") && min(y_i) <= 0 && !delta) {
      cli_abort("Gamma and lognormal must have response values > 0.")
    }
  }
  # For distribution_column cases, defer validation until after processing
}

#' Process binomial response variable and size for single family case
#'
#' @param family Single family object
#' @param mf Model frame list
#' @param X_ij Design matrix list
#' @param delta Logical indicating delta model
#' @param distribution_column Character string or NULL
#' @param weights Optional weights vector
#' @return List with processed y_i and size vectors
#' @keywords internal
#' @noRd
process_binomial_response <- function(family, mf, X_ij, delta, distribution_column, weights = NULL) {
  family_info <- safe_family_check(family, distribution_column)
  
  # Initialize default values
  size <- rep(1, nrow(X_ij[[1]]))
  y_i <- NULL
  
  if (family_info$is_single && identical(family_info$first_family, "binomial") && !delta) {
    ## call this to catch the factor / matrix cases
    y_i <- model.response(mf[[1]], type = "any")
    ## allow character
    if (is.character(y_i)) {
      y_i <- model.response(mf[[1]], type = "factor")
      if (nlevels(y_i) > 2) {
        cli_abort("More than 2 levels detected for response")
      }
    }
    if (is.factor(y_i)) {
      if (nlevels(y_i) > 2) {
        cli_abort("More than 2 levels detected for response")
      }
      ## following glm, 'success' is interpreted as the factor not
      ## having the first level (and hence usually of having the
      ## second level).
      y_i <- pmin(as.numeric(y_i) - 1, 1)
      size <- rep(1, length(y_i))
    } else {
      if (is.matrix(y_i)) { # yobs=cbind(success, failure)
        size <- y_i[, 1] + y_i[, 2]
        yobs <- y_i[, 1] # successes
        y_i <- yobs
      } else {
        if (all(y_i %in% c(0, 1))) { # binary
          size <- rep(1, length(y_i))
        } else { # proportions
          if (!is.null(weights)) {
            y_i <- weights * y_i
            size <- weights
            weights <- rep(1, length(y_i))
          } else {
            size <- rep(1, length(y_i))
          }
        }
      }
    }
    
    # https://github.com/pbs-assess/sdmTMB/issues/172
    if (is.logical(y_i)) {
      msg <- paste0(
        "We recommend against using `TRUE`/`FALSE` ",
        "response values if you are going to use the `visreg::visreg()` ",
        "function after. Consider converting to integer with `as.integer()`."
      )
      cli_warn(msg)
    }
  }
  
  return(list(y_i = y_i, size = size, weights = weights))
}

#' Validate late family constraints after distribution processing
#'
#' @param family_info Safe family info from safe_family_check()
#' @param dist_result Result from process_distribution_column()
#' @param delta Logical indicating delta model
#' @keywords internal
#' @noRd
validate_late_family_constraints <- function(family_info, dist_result = NULL, delta) {
  if (family_info$is_single) {
    family_obj <- family_info$family_obj
    first_fam <- family_info$first_family
    
    # Student df parameter
    df <- if (first_fam == "student" && "df" %in% names(family_obj)) family_obj$df else 3
    
    # Return values needed by fit.R
    return(list(
      df = df,
      first_family = first_fam,
      family_obj = family_obj
    ))
  } else {
    # For distribution_column cases, use processed results
    return(list(
      df = 3,  # Default df for mixed cases
      first_family = NULL,
      family_obj = NULL
    ))
  }
}

#' Transform response matrix for mixed delta/non-delta families
#'
#' @param y_i Original response vector
#' @param delta Logical indicating if single family is delta
#' @param distribution_column Character string of column name, or NULL
#' @param dist_result Result from process_distribution_column()
#' @return Transformed response matrix
#' @keywords internal
#' @noRd
transform_response_matrix <- function(y_i, delta, distribution_column, dist_result) {
  if (is.null(distribution_column)) {
    # Single family case - use original logic
    if (delta) {
      return(cbind(ifelse(y_i > 0, 1, 0), ifelse(y_i > 0, y_i, NA_real_)))
    } else {
      return(matrix(y_i, ncol = 1L))
    }
  } else {
    # Mixed family case - vectorized transformation based on component usage
    n_obs <- length(y_i)
    n_m_effective <- dist_result$n_m_effective
    
    # Create logical vector for delta vs non-delta observations
    is_delta_obs <- dist_result$components_per_obs == 2L
    
    # Initialize transformed response matrix
    y_i_transformed <- array(NA_real_, dim = c(n_obs, n_m_effective))
    
    # Vectorized assignment for first component
    y_i_transformed[, 1] <- ifelse(is_delta_obs, 
                                   ifelse(y_i > 0, 1, 0),  # Delta: 0/1 indicator
                                   y_i)                     # Non-delta: original values
    
    # For second component (if needed)
    if (n_m_effective > 1L) {
      y_i_transformed[, 2] <- ifelse(is_delta_obs,
                                     ifelse(y_i > 0, y_i, NA_real_),  # Delta: positive values
                                     NA_real_)                        # Non-delta: unused
    }
    
    return(y_i_transformed)
  }
}
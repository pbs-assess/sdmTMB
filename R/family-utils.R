#' Build Family Arrays for TMB
#'
#' Internal function to create vectorized family/link arrays for both single-family
#' and integrated (distribution_column) models. This function encapsulates all
#' family-related processing in one place.
#'
#' @param family Either a single family object (e.g., tweedie()) or a named list
#'   of family objects for integrated models
#' @param distribution_column Character string of column name in data specifying
#'   distribution per observation, or NULL for single-family models
#' @param data Data frame containing the distribution column (if applicable)
#' @param n_obs Number of observations
#'
#' @return List containing:
#'   \item{family}{Integer array [n_obs x n_m] of family codes}
#'   \item{link}{Integer array [n_obs x n_m] of link codes}
#'   \item{e_i}{Integer vector [n_obs] mapping observations to phi parameter groups (0-indexed)}
#'   \item{poisson_link_delta}{Integer vector [n_obs] indicating poisson-link-delta per observation}
#'   \item{n_m_effective}{Integer, effective number of model components (1 or 2)}
#'   \item{component_usage}{Integer array [n_obs x n_m] indicating which components are active}
#'   \item{n_families}{Integer, number of unique families}
#'   \item{family_list}{List of family objects (for downstream use)}
#'   \item{components}{Integer vector [n_families] indicating components per family}
#'
#' @keywords internal
#' @noRd
build_family_arrays <- function(family, distribution_column, data, n_obs) {

  # ---- Determine single vs integrated case ----
  if (is.null(distribution_column)) {
    # Single family case
    if (is.list(family) && !is.null(family$family)) {
      # Standard family object
      n_families <- 1L
      family_list <- list(family)
      family_indices <- rep(1L, n_obs)
    } else {
      cli_abort("Invalid `family` argument. Expected a family object like tweedie() or gaussian().")
    }
  } else {
    # Integrated model case - match tinyVAST interface
    if (!is.list(family) || is.null(names(family))) {
      cli_abort("When `distribution_column` is provided, `family` must be a named list of family objects.")
    }

    # Extract distribution values from data
    if (!distribution_column %in% names(data)) {
      cli_abort("Column '{distribution_column}' not found in data.")
    }
    dist_values <- data[[distribution_column]]

    # Match observations to family names - following tinyVAST exactly
    family_indices <- match(dist_values, names(family))

    # Validate - match tinyVAST error message
    if (any(is.na(family_indices))) {
      missing <- unique(dist_values[is.na(family_indices)])
      cli_abort("`data[,distribution_column]` has values that don't match `names(family)`: {paste(missing, collapse = ', ')}")
    }

    n_families <- length(family)
    family_list <- family
  }

  # ---- Determine components per family ----
  # Following tinyVAST: delta families have 2 components, others have 1
  components <- sapply(family_list, function(x) {
    if (isTRUE(x$delta)) 2L else 1L
  })

  # Maximum components needed determines n_m
  n_m_effective <- max(components)

  # ---- Build e_i mapping for phi parameters (0-indexed for TMB) ----
  # Maps each observation to its family group
  e_i <- as.integer(family_indices) - 1L

  # ---- Initialize output arrays ----
  family_array <- array(NA_integer_, dim = c(n_obs, n_m_effective))
  link_array <- array(NA_integer_, dim = c(n_obs, n_m_effective))
  component_usage <- array(0L, dim = c(n_obs, n_m_effective))
  poisson_link_delta <- integer(n_obs)

  # ---- Fill arrays by observation ----
  # Vectorized by observation
  for (i in seq_len(n_obs)) {
    fam_idx <- family_indices[i]
    fam_obj <- family_list[[fam_idx]]
    n_comp <- components[fam_idx]

    # Fill family and link codes for active components
    for (m in seq_len(n_comp)) {
      family_array[i, m] <- .valid_family[fam_obj$family[m]]
      link_array[i, m] <- .valid_link[fam_obj$link[m]]
      component_usage[i, m] <- 1L
    }

    # Check for poisson-link delta
    poisson_link_delta[i] <- as.integer(isTRUE(fam_obj$type == "poisson_link_delta"))
  }

  # ---- Return structured list ----
  return(list(
    family = family_array,
    link = link_array,
    e_i = e_i,
    poisson_link_delta = poisson_link_delta,
    n_m_effective = n_m_effective,
    component_usage = component_usage,
    n_families = n_families,
    family_list = family_list,
    components = components
  ))
}

#' Turn sdmTMB model output into a tidy data frame
#'
#' @param x Output from [sdmTMB()].
#' @param effects A character value. One of `"fixed"` ('fixed' or main-effect
#'   parameters), `"ran_pars"` (standard deviations, spatial range, and other
#'   random effect and dispersion-related terms), `"ran_vals"` (individual
#'   random intercepts or slopes, if included; behaves like `ranef()`), or `"ran_vcov"` (list
#'   of variance covariance matrices for the random effects, by model and group).
#' @param conf.int Include a confidence interval?
#' @param conf.level Confidence level for CI.
#' @param exponentiate Whether to exponentiate the fixed-effect coefficient
#'   estimates and confidence intervals.
#' @param model Which model to tidy if a delta model (1 or 2). The `model` will be
#'   ignored when effects is `"ran_vals"` (all returned in a single dataframe)
#'
#' @param silent Omit any messages?
#' @param ... Extra arguments (not used).
#'
#' @return A data frame
#' @details
#' Follows the conventions of the \pkg{broom} and \pkg{broom.mixed} packages.
#'
#' Currently, `effects = "ran_pars"` also includes dispersion-related terms
#' (e.g., `phi`), which are not actually associated with random effects.
#'
#' Standard errors for spatial variance terms fit in log space (e.g., variance
#' terms, range, or parameters associated with the observation error) are
#' omitted to avoid confusion. Confidence intervals are still available.
#'
#' @export
#'
#' @importFrom assertthat assert_that
#' @importFrom stats plogis
#' @examples
#' fit <- sdmTMB(density ~ poly(depth_scaled, 2, raw = TRUE),
#'   data = pcod_2011, mesh = pcod_mesh_2011,
#'   family = tweedie()
#' )
#' tidy(fit)
#' tidy(fit, conf.int = TRUE)
#' tidy(fit, "ran_pars", conf.int = TRUE)
#'
#' pcod_2011$fyear <- as.factor(pcod_2011$year)
#' fit <- sdmTMB(density ~ poly(depth_scaled, 2, raw = TRUE) + (1 | fyear),
#'   data = pcod_2011, mesh = pcod_mesh_2011,
#'   family = tweedie()
#' )
#' tidy(fit, "ran_vals")

tidy.sdmTMB <- function(x, effects = c("fixed", "ran_pars", "ran_vals", "ran_vcov"), model = 1,
                 conf.int = TRUE, conf.level = 0.95, exponentiate = FALSE,
                 silent = FALSE, ...) {
  effects <- match.arg(effects)
  assert_that(is.logical(exponentiate))
  assert_that(is.logical(conf.int))
  if (conf.int) {
    assert_that(is.numeric(conf.level),
      conf.level > 0, conf.level < 1,
      length(conf.level) == 1,
      msg = "`conf.level` must be length 1 and between 0 and 1")
  }

  crit <- stats::qnorm(1 - (1 - conf.level) / 2)
  if (exponentiate) trans <- exp else trans <- I

  reinitialize(x)

  delta <- isTRUE(x$family$delta)
  assert_that(is.numeric(model))
  assert_that(length(model) == 1L)
  if (delta) assert_that(model %in% c(1, 2), msg = "`model` must be 1 or 2.")
  if (!delta) assert_that(model == 1, msg = "Only one model: `model` must be 1.")

  se_rep <- as.list(x$sd_report, "Std. Error", report = TRUE)
  est_rep <- as.list(x$sd_report, "Estimate", report = TRUE)
  se <- as.list(x$sd_report, "Std. Error", report = FALSE)
  est <- as.list(x$sd_report, "Estimate", report = FALSE)

  se <- c(se, se_rep)
  est <- c(est, est_rep)
  # cleanup:
  est$epsilon_st <- NULL
  est$zeta_s <- NULL
  est$omega_s <- NULL
  est$ln_H_input <- NULL

  se$epsilon_st <- NULL
  se$zeta_s <- NULL
  se$omega_s <- NULL
  se$ln_H_input <- NULL

  subset_pars <- function(p, model) {
    p$b_j <- if (model == 1) p$b_j else p$b_j2
    p$ln_tau_O <- p$ln_tau_O[model]
    p$ln_tau_Z <- p$ln_tau_Z[model]
    p$ln_tau_E <- p$ln_tau_E[model]
    p$ln_kappa <- as.numeric(p$ln_kappa[,model])
    p$ln_phi <- p$ln_phi[model]
    p$ln_tau_V <- as.numeric(p$ln_tau_V[,model])
    p$ar1_phi <- as.numeric(p$ar1_phi[model])
    p$log_sigma_O <- as.numeric(p$log_sigma_O[1,model])
    p$log_sigma_E <- as.numeric(p$log_sigma_E[1,model])
    p$log_sigma_Z <- as.numeric(p$log_sigma_Z[,model])
    p$log_range <- as.numeric(p$log_range[,model])

    p$phi <- p$phi[model]
    p$range <- as.numeric(p$range[,model])
    p$sigma_E <- as.numeric(p$sigma_E[1,model])
    p$sigma_O <- as.numeric(p$sigma_O[1,model])
    p$sigma_Z <- as.numeric(p$sigma_Z[,model])

    # if delta, a single AR1 -> rho_time_unscaled is a 1x2 matrix
    p$rho_time <- 2 * plogis(p$rho_time_unscaled[,model]) - 1
    p
  }

  est <- subset_pars(est, model)
  se <- subset_pars(se, model)

  if (x$family$family[[model]] %in% c("binomial", "poisson")) {
    se$ln_phi <- NULL
    est$ln_phi <- NULL
    se$phi <- NULL
    est$phi <- NULL
  }

  ii <- 1

  # grab fixed effects:
  .formula <- x$split_formula[[model]]$form_no_bars
  .formula <- remove_s_and_t2(.formula)
  if (!"mgcv" %in% names(x)) x[["mgcv"]] <- FALSE
  fe_names <- colnames(model.matrix(.formula, x$data))

  b_j <- est$b_j[!fe_names == "offset", drop = TRUE]
  b_j_se <- se$b_j[!fe_names == "offset", drop = TRUE]
  fe_names <- fe_names[!fe_names == "offset"]
  out <- data.frame(term = fe_names, estimate = b_j, std.error = b_j_se, stringsAsFactors = FALSE)

  if (x$tmb_data$threshold_func > 0) {
    if (x$threshold_function == 1L) {
      par_name <- paste0(x$threshold_parameter, c("-slope", "-breakpt"))
      estimates <- est$b_threshold[,model,drop=TRUE]
      ses <- se$b_threshold[,model,drop=TRUE]
    } else {
      par_name <- paste0(x$threshold_parameter, c("-s50", "-s95", "-smax"))
      estimates <- c(est$s50[model], est$s95[model], est$s_max[model])
      ses <- c(se$s50[model], se$s95[model], se$s_max[model])
    }
    out <- rbind(
      out,
      data.frame(
        term = par_name, estimate = estimates,
        std.error = ses, stringsAsFactors = FALSE
      )
    )
  }

  if (x$tmb_data$has_smooths) {
    p <- print_smooth_effects(x, silent = FALSE, m = model)
    mm <- p$smooth_effects
    out <- rbind(
      out,
      data.frame(
        term = rownames(mm),
        estimate = mm[, "bs"],
        std.error = mm[, "bs_se"],
        stringsAsFactors = FALSE
      )
    )
  }

  if (conf.int) {
    out$conf.low <- as.numeric(trans(out$estimate - crit * out$std.error))
    out$conf.high <- as.numeric(trans(out$estimate + crit * out$std.error))
  }
  # must wrap in as.numeric() otherwise I() leaves 'AsIs' class that affects emmeans package
  out$estimate <- as.numeric(trans(out$estimate))
  if (exponentiate) out$std.error <- NULL

  out_re <- list()
  log_name <- c("log_range")
  name <- c("range")
  if (!isTRUE(is.na(x$tmb_map$ln_phi))) {
    log_name <- c(log_name, "ln_phi")
    name <- c(name, "phi")
  }
  if (x$tmb_data$include_spatial[model]) {
    log_name <- c(log_name, "log_sigma_O")
    name <- c(name, "sigma_O")
  }
  if (!x$tmb_data$spatial_only[model]) {
    log_name <- c(log_name, "log_sigma_E")
    name <- c(name, "sigma_E")
  }
  if (x$tmb_data$spatial_covariate) {
    log_name <- c(log_name, "log_sigma_Z")
    name <- c(name, "sigma_Z")
  }
  if (x$tmb_data$random_walk) {
    log_name <- c(log_name, "ln_tau_V")
    name <- c(name, "tau_V")
  }
  if (!all(est$rho_time == 0)) {
    log_name <- c(log_name, "rho_time_unscaled")
    name <- c(name, "rho_time")
  }

  j <- 0
  if (!"log_range" %in% names(est)) {
    cli_warn("This model was fit with an old version of sdmTMB. Some parameters may not be available to the tidy() method. Re-fit the model with the current version of sdmTMB if you need access to any missing parameters.")
  }

  for (i in name) {
    j <- j + 1
    if (i %in% names(est)) {
      .e <- est[[log_name[j]]]
      .se <- se[[log_name[j]]]
      .e <- if (is.null(.e)) NA else .e
      .se <- if (is.null(.se)) NA else .se

      non_log_name <- gsub("ln_", "", gsub("log_", "", log_name))
      this <- non_log_name[j]
      if (this == "tau_V") this <- "sigma_V"
      if (this == "rho_time_unscaled") this <- "rho_time"

      this_se <- as.numeric(se[[this]])
      this_est <- as.numeric(est[[this]])
      if (length(this_est) && !(all(this_se == 0) && all(this_est == 0))) {
        out_re[[i]] <- data.frame(
          term = i, estimate = this_est, std.error = this_se,
          conf.low = exp(.e - crit * .se),
          conf.high = exp(.e + crit * .se),
          stringsAsFactors = FALSE
        )
        if(this == "rho_time") {
          out_re[[i]] <- data.frame(
            term = i,
            estimate = this_est,
            # use delta method to get SE in normal space
            std.error = 2 * plogis (.e[,model]) * (1 - plogis (.e[,model])) * .se[,model],
            # don't use delta-method for CIs, because they can be outside (-1,1)
            conf.low = 2 * plogis(.e[,model] - crit * .se[,model]) - 1,
            conf.high = 2 * plogis(.e[,model] + crit * .se[,model]) - 1,
            stringsAsFactors = FALSE
          )

        }
      }
      ii <- ii + 1
    }
  }
  discard <- unlist(lapply(out_re, function(x) length(x) == 1L)) # e.g. old models and phi
  out_re[discard] <- NULL

  if ("tweedie" %in% x$family$family) {
    out_re$tweedie_p <- data.frame(
      term = "tweedie_p", estimate = plogis(est$thetaf) + 1,
      std.error = se$tweedie_p, stringsAsFactors = FALSE)
    out_re$tweedie_p$conf.low <- plogis(est$thetaf - crit * se$thetaf) + 1
    out_re$tweedie_p$conf.high <- plogis(est$thetaf + crit * se$thetaf) + 1
    ii <- ii + 1
  }

  if ("ar1_phi" %in% names(est)) {
    ar_phi <- est$ar1_phi
    ar_phi_se <- se$ar1_phi
    rho_est <- 2 * stats::plogis(ar_phi) - 1
    rho_lwr <- 2 * stats::plogis(ar_phi - crit * ar_phi_se) - 1
    rho_upr <- 2 * stats::plogis(ar_phi + crit * ar_phi_se) - 1
    out_re[[ii]] <- data.frame(
      term = "rho", estimate = rho_est, std.error = NA,
      conf.low = rho_lwr, conf.high = rho_upr, stringsAsFactors = FALSE
    )
    ii <- ii + 1
  }

  if (all(!x$tmb_data$include_spatial) && all(x$tmb_data$spatial_only)) out_re$range <- NULL

  out_re <- do.call("rbind", out_re)
  row.names(out_re) <- NULL

  if (identical(est$ln_tau_E, 0)) out_re <- out_re[out_re$term != "sigma_E", ]
  if (identical(est$ln_tau_V, 0)) out_re <- out_re[out_re$term != "sigma_V", ]
  if (identical(est$ln_tau_O, 0)) out_re <- out_re[out_re$term != "sigma_O", ]
  if (identical(est$ln_tau_Z, 0)) out_re <- out_re[out_re$term != "sigma_Z", ]
  if (is.na(x$tmb_map$ar1_phi[model])) out_re <- out_re[out_re$term != "rho", ]

  if (!conf.int) {
    out_re[["conf.low"]] <- NULL
    out_re[["conf.high"]] <- NULL
  }

  if (sum(x$tmb_data$n_re_groups) > 0L) { # we have random intercepts/slopes
    temp <- get_re_tidy_list(x, crit = crit, model = model)
    cov_mat_list <- list(est = temp$cov_matrices)
    if(conf.int) {
      cov_mat_list[["lo"]] <- temp$cov_matrices_lo
      cov_mat_list[["hi"]] <- temp$cov_matrices_hi
    }
    out_ranef <- temp$out_ranef
  } else {
    cov_mat_list <- NULL
    out_ranef <- NULL
  }

  # optional time-varying random components
  if (!is.null(x$time_varying)) {
    tv_names <- colnames(model.matrix(x$time_varying, x$data))
    time_slices <- x$time_lu$time_from_data
    yrs <- rep(time_slices, times = length(tv_names))

    out_ranef_tv <- data.frame(
      model = model,
      term = paste0(rep(tv_names, each = length(time_slices)), ":", yrs),
      estimate = c(est$b_rw_t),
      std.error = c(se$b_rw_t),
      conf.low = c(est$b_rw_t) - crit * c(se$b_rw_t),
      conf.high = c(est$b_rw_t) + crit * c(se$b_rw_t),
      stringsAsFactors = FALSE
    )

    if(is.null(out_ranef)) {
      out_ranef <- out_ranef_tv
    } else {
      out_ranef_tv$group_name <- NA
      out_ranef_tv$level_ids <- NA
      out_ranef <- rbind(out_ranef, out_ranef_tv)
    }
  }

  out <- unique(out) # range can be duplicated
  out_re <- unique(out_re)

  if (requireNamespace("tibble", quietly = TRUE)) {
    frm <- tibble::as_tibble
  } else {
    frm <- as.data.frame
  }

  if (effects == "fixed") {
    return(frm(out))
  } else if (effects == "ran_vals") {
    return(frm(out_ranef))
  } else if (effects == "ran_pars") {
    return(frm(out_re))
  } else if (effects == "ran_vcov") {
    return(cov_mat_list)
  } else {
    cli_abort("The specified 'effects' type is not available.")
  }
}

# Extract and format random effect estimates from sdmTMB model output
#
# This function extracts random effect estimates, including individual random intercepts
# and slopes, as well as covariance matrices, from an `sdmTMB` model output. It formats
# them into a structured list for further analysis.
#
# @param x An `sdmTMB` model object containing estimated random effects.
# @param crit The critical value for confidence interval computation,
#   typically derived from a normal distribution quantile, e.g.
#   crit = stats::qnorm(1 - (1 - conf.level) / 2)
# @importFrom stats aggregate

get_re_tidy_list <- function(x, crit, model = 1) {
  re_b_dfs <- add_model_index(x$split_formula, "re_b_df")
  re_b_df <- do.call(rbind, re_b_dfs)
  names(re_b_df)[names(re_b_df) == "group_indices"] <- "group_id"

  # this function just expands each row from start: end
  expand_row <- function(level_id, start, end, group_id, model) {
    seq_len <- end - start + 1
    data.frame(
      level_ids = rep(level_id, seq_len),
      index = seq(from = start, to = end),
      group_id = rep(group_id, seq_len),
      model = rep(model, seq_len)
    )
  }

  # apply to each row and combine the results
  expanded_rows <- Map(
    expand_row, re_b_df$level_ids,
    re_b_df$start, re_b_df$end,
    re_b_df$group_id, re_b_df$model
  )

  re_b_df <- do.call(rbind, expanded_rows) # list to df
  rownames(re_b_df) <- NULL # reset row names

  # this is all as before
  re_indx <- grep("re_b_pars", names(x$sd_report$value), fixed = TRUE)
  non_nas <- !is.na(x$tmb_map$re_b_pars) # parameters that don't get mapped off

  re_b_df$estimate <- x$sd_report$value[re_indx][non_nas]
  re_b_df$std.error <- x$sd_report$sd[re_indx][non_nas]
  re_b_df$conf.low <- re_b_df$estimate - crit * re_b_df$std.error
  re_b_df$conf.high <- re_b_df$estimate + crit * re_b_df$std.error
  re_b_df$index <- NULL
  re_b_df$group_name <- NA
  re_b_df$term <- NA
  for (i in seq_len(length(x$split_formula))) {
    groupnames <- names(x$split_formula[[i]]$re_cov_terms$cnms)
    for (j in seq_len(length(x$split_formula[[i]]$barnames))) {
      model_grp <- which(re_b_df$model == i & re_b_df$group_id == j)
      re_b_df$group_name[model_grp] <- groupnames[j]
      re_b_df$term[model_grp] <- rep(x$split_formula[[i]]$re_cov_terms$cnms[[j]], length.out = length(model_grp))
    }
  }
  group_key <- stats::aggregate(group_name ~ model + group_id, data = re_b_df, FUN = function(x) x[1])

  # more sensible re-ordering
  re_b_df$group_id <- NULL
  re_b_df <- re_b_df[, c("model", "group_name", "term", "level_ids", "estimate", "std.error", "conf.low", "conf.high")]
  # remove ":" in the level_ids
  re_b_df$level_ids <- sapply(strsplit(re_b_df$level_ids, ":"), function(x) x[2])
  out_ranef <- re_b_df
  row.names(out_ranef) <- NULL

  re_cov_dfs <- add_model_index(x$split_formula, "re_df")
  # add group names
  re_cov_dfs <- lapply(re_cov_dfs, function(df) {
    df$group <- row.names(df)
    df
  })
  re_cov_df <- do.call(rbind, re_cov_dfs)
  re_cov_df$rows <- re_cov_df$rows + 1 # increement rows/cols from 1, not 0
  re_cov_df$cols <- re_cov_df$cols + 1
  row.names(re_cov_df) <- NULL

  # make sure group index is included correctly
  for (i in seq_len(nrow(re_cov_df))) {
    indx <- which(group_key$model == re_cov_df$model[i] & group_key$group_id == (re_cov_df$group_indices[i] + 1))
    re_cov_df$group[i] <- group_key$group_name[indx]
  }

  re_indx <- grep("re_cov_pars", names(x$sd_report$value), fixed = TRUE)
  non_nas <- which(x$sd_report$value[re_indx] != 0) # remove parameter that gets mapped off
  re_cov_df$estimate <- x$sd_report$value[re_indx][non_nas]
  re_cov_df$std.error <- x$sd_report$sd[re_indx][non_nas]
  re_cov_df$conf.low <- re_cov_df$estimate - crit * re_cov_df$std.error
  re_cov_df$conf.high <- re_cov_df$estimate + crit * re_cov_df$std.error
  # the SD parameters are returned in log space -- use delta method to generate CIs
  est <- exp(re_cov_df$estimate) # estimate in normal space
  sd_est <- est * re_cov_df$std.error # SE in normal space
  sds <- which(re_cov_df$is_sd == 1) # index which elements are SDs
  re_cov_df$estimate[sds] <- est[sds]
  re_cov_df$conf.low[sds] <- est[sds] - crit * sd_est[sds]
  re_cov_df$conf.high[sds] <- est[sds] + crit * sd_est[sds]

  re_cov_df <- re_cov_df[, c("rows", "cols", "model", "group", "estimate", "std.error", "conf.low", "conf.high")]
  cov_matrices_lo = create_cov_matrices(re_cov_df, col_name = "conf.low", model = model)
  cov_matrices_hi = create_cov_matrices(re_cov_df, col_name = "conf.high", model = model)
  cov_matrices_est = create_cov_matrices(re_cov_df, model = model)
  list(
    out_ranef = out_ranef,
    cov_matrices = cov_matrices_est,
    cov_matrices_lo = cov_matrices_lo,
    cov_matrices_hi = cov_matrices_hi
  )
}

create_cov_matrices <- function(df, col_name = "estimate", model = 1) {
  # Initialize an empty list to store the covariance matrices
  cov_matrices <- list()

  # Group by model and group, and then create a covariance matrix for each group
  for (m in model) {
    model_df <- df[df$model == m,,drop=FALSE]
    for (group in unique(model_df$group)) {
      group_df <- model_df[model_df$group == group,,drop=FALSE]
      # Create an empty matrix
      cov_matrix <- matrix(NA_real_, nrow = max(group_df$rows), ncol = max(group_df$cols))
      # Fill the matrix with the estimates
      for (i in seq_len(nrow(group_df))) {
        cov_matrix[group_df$rows[i], group_df$cols[i]] <- group_df[[col_name]][i]
      }
      # Add the matrix to the list
      cov_matrices[[paste("Model", m, "Group", group)]] <- cov_matrix
    }
  }
  cov_matrices
}

#' @importFrom generics tidy
#' @export
generics::tidy

#' @rdname tidy.sdmTMB
#' @export
tidy.sdmTMB_cv <- function(x, ...) {
  x <- x$models
  out <- lapply(seq_along(x), function(i) {
    df <- tidy.sdmTMB(x[[i]], ...)
    df$cv_split <- i # add a model index column
    df
  })
  do.call("rbind", out)
}

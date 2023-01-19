# from glmmTMB:
named_list <- function(...) {
  L <- list(...)
  snm <- sapply(substitute(list(...)), deparse)[-1L]
  if (is.null(nm <- names(L))) {
    nm <- snm
  }
  if (any(nonames <- nm == "")) {
    nm[nonames] <- snm[nonames]
  }
  stats::setNames(L, nm)
}

print_model_info <- function(x) {
  delta <- isTRUE(x$family$delta)
  spatial_only <- as.logical(x$tmb_data$spatial_only)
  fit_by <- if (isTRUE(x$reml)) "REML" else "ML"
  if (all(spatial_only)) {
    title <- paste0("Spatial model fit by ", fit_by, " ['sdmTMB']\n")
  } else {
    title <- paste0("Spatiotemporal model fit by ", fit_by, " ['sdmTMB']\n")
  }
  if (all(x$spatial == "off") && all(x$spatiotemporal == "off")) {
    title <- paste0("Model fit by ", fit_by, " ['sdmTMB']\n")
  }

  aniso <- as.logical(x$tmb_data$anisotropy)
  # covariance <- paste0("Covariance: ", if (aniso) "anisotropic" else "isotropic", "\n")
  covariance <- paste0(if (aniso) "anisotropic" else "isotropic")

  formula <- paste0("Formula: ", deparse(x$call$formula), "\n")

  if (deparse(x$call$time) != "NULL") {
    time <- paste0("Time column: ", deparse(x$call$time), "\n")
    time <- gsub('\\"', "", time)
    time <- gsub("\\'", "", time)
  } else {
    time <- NULL
  }

  mesh <- paste0("Mesh: ", deparse(x$call$mesh), " (", covariance, " covariance)\n")
  data <- paste0("Data: ", deparse(x$call$data), "\n")

  # From sdmTMB_cv():
  if (length(mesh) > 1L) mesh <- NULL
  if (length(data) > 1L) data <- NULL

  if ("clean_name" %in% names(x$family)) {
    overall_family <- x$family$clean_name
  } else {
    overall_family <- paste0(x$family$family[1], "(link = '", x$family$link[1], "')")
  }
  overall_family <- paste0("Family: ", overall_family, "\n")

  if (delta) {
    family1 <- paste0("Family: ", x$family$family[1], "(link = '", x$family$link[1], "')")
    family2 <- paste0("Family: ", x$family$family[2], "(link = '", x$family$link[2], "')")
  } else {
    family1 <- NULL
    family2 <- NULL
  }

  criterion <- paste0(fit_by, " criterion at convergence: ", mround(x$model$objective, 3), "\n")

  named_list(
    delta,
    spatial_only,
    title,
    time,
    mesh,
    formula,
    data,
    family1,
    family2,
    overall_family,
    criterion,
    covariance
  )
}

print_main_effects <- function(x, m = 1) {
  b <- as.data.frame(tidy(x, model = m, silent = TRUE))
  b$estimate <- round(b$estimate, 2L)
  b$std.error <- round(b$std.error, 2L)
  mm <- cbind(b$estimate, b$std.error)
  colnames(mm) <- c("coef.est", "coef.se")
  row.names(mm) <- b$term
  mm
}

print_smooth_effects <- function(x, m = 1) {
  sr <- x$sd_report
  sr_se <- as.list(sr, "Std. Error")
  sr_est <- as.list(sr, "Estimate")
  if (x$tmb_data$has_smooths) {
    bs_se <- round(sr_se$bs[, m], 2L)
    bs <- round(sr_est$bs[, m], 2L)
    sm <- parse_smoothers(formula = x$formula[[m]], data = x$data)
    sm_names <- unlist(lapply(sm$Zs, function(x) attr(x, "s.label")))
    sm_names <- gsub("\\)$", "", gsub("s\\(", "", sm_names))
    sm_names_bs <- paste0("s", sm_names)
    sm_classes <- unlist(sm$classes)
    xx <- lapply(sm_names_bs, function(.x) { # split out 2D + smooths
      n_sm <- grep(",", .x) + 1
      if (length(n_sm)) {
        .x <- paste(.x, seq_len(n_sm), sep = "_")
        gsub(",", "", .x)
      } else {
        .x
      }
    })
    sm_names_bs <- unlist(xx)
    sm_names_sds <- paste0("sds(", sm_names, ")")
    mm_sm <- cbind(bs, bs_se)

    .sm_names_bs <- sm_names_bs[!sm_classes %in% c("cc.smooth.spec", "cs.smooth.spec")]
    if (length(.sm_names_bs) != nrow(mm_sm)) {
      if (nrow(mm_sm) > 0L) {
        msg <- c(
          "Smoother fixed effect matrix names could not be retrieved.",
          "This could be from setting m != 2 in s() or other unanticipated s() arguments.",
          "This does not affect model fitting.",
          "We'll use generic covariate names ('scovariate') here intead."
        )
        cli_warn(msg)
        row.names(mm_sm) <- paste0("scovariate-", seq_len(nrow(mm_sm)))
      } else {
        mm_sm <- NULL
      }
    } else {
      row.names(mm_sm) <- .sm_names_bs
    }

    smooth_sds <- round(exp(sr_est$ln_smooth_sigma[, m]), 2L)
    re_sm_mat <- matrix(NA_real_, nrow = length(smooth_sds), ncol = 1L)
    re_sm_mat[, 1] <- smooth_sds
    rownames(re_sm_mat) <- sm_names_sds
    colnames(re_sm_mat) <- "Std. Dev."
  } else {
    re_sm_mat <- NULL
    mm_sm <- NULL
  }
  list(smooth_effects = mm_sm, smooth_sds = re_sm_mat)
}

print_iid_re <- function(x, m = 1) {
  .tidy <- tidy(x, "ran_pars", model = m, silent = TRUE)
  if ("sigma_G" %in% .tidy$term) {
    re_int_names <- barnames(x$split_formula[[1]]$reTrmFormulas) # TODO DELTA HARDCODED TO 1
    re_int_mat <- matrix(NA_real_, nrow = length(re_int_names), ncol = 1L)
    re_int_mat[, 1L] <- round(.tidy$estimate[.tidy$term == "sigma_G"], 2L)
    rownames(re_int_mat) <- re_int_names
    colnames(re_int_mat) <- "Std. Dev."
  } else {
    re_int_mat <- NULL
  }
  re_int_mat
}

print_time_varying <- function(x, m = 1) {
  if (!is.null(x$time_varying)) {
    sr <- x$sd_report
    sr_se <- as.list(sr, "Std. Error")
    sr_est <- as.list(sr, "Estimate")
    b_rw_t_est <- sr_est$b_rw_t[, , m]
    b_rw_t_se <- sr_se$b_rw_t[, , m]
    tv_names <- colnames(model.matrix(x$time_varying, x$data))
    mm_tv <- cbind(round(as.numeric(b_rw_t_est), 2L), round(as.numeric(b_rw_t_se), 2L))
    colnames(mm_tv) <- c("coef.est", "coef.se")
    time_slices <- sort(unique(x$data[[x$time]]))
    row.names(mm_tv) <- paste(rep(tv_names, each = length(time_slices)), time_slices, sep = "-")
  } else {
    mm_tv <- NULL
  }
  mm_tv
}

print_range <- function(x, m = 1) {
  b <- tidy(x, effects = "ran_pars", model = m, silent = TRUE)
  range <- b$estimate[b$term == "range"]
  if (is.null(range)) {
    return(NULL)
  }

  # if (as.logical(x$tmb_data$anisotropy)) {
  #   aa <- plot_anisotropy(x, return_data = TRUE)
  #   # range <- aa$a[aa$random_field == "spatial"]
  # }

  range <- mround(range, 2L)
  range_text <- if (x$tmb_data$share_range[m]) {
    paste0("Matern range: ", range[1], "\n")
  } else {
    paste0(
      "Matern range (spatial): ", range[1], "\n",
      "Matern range (spatiotemporal): ", range[2], "\n"
    )
  }
  if (x$spatial[m] == "off" && x$spatiotemporal[m] == "off") {
    range_text <- NULL
  }

  if (as.logical(x$tmb_data$anisotropy)) {
    range_text <- "Matern range: anisotropic covariance; see plot_anisotropy().\n"
  }
  range_text
}

print_other_parameters <- function(x, m = 1) {
  b <- tidy(x, "ran_pars", model = m, silent = TRUE)

  get_term_text <- function(term_name = "", pretext = "") {
    if (term_name %in% b$term) {
      a <- mround(b$estimate[b$term == term_name], 2L)
      a <- paste0(pretext, ": ", a, "\n")
    } else {
      a <- ""
    }
    a
  }

  phi <- get_term_text("phi", "Dispersion parameter")
  tweedie_p <- get_term_text("tweedie_p", "Tweedie p")
  sigma_O <- get_term_text("sigma_O", "Spatial SD")
  sigma_E <- get_term_text("sigma_E", "Spatiotemporal SD")
  rho <- get_term_text("rho", "Spatiotemporal AR1 correlation (rho)")

  if ("sigma_Z" %in% b$term) {
    a <- mround(b$estimate[b$term == "sigma_Z"], 2L)
    sigma_Z <- paste0("Spatially varying coefficient SD (", x$spatial_varying,  "): ", a, "\n")
    sigma_Z <- paste(sigma_Z, collapse = "")
  } else {
    sigma_Z <- ""
  }

  named_list(phi, tweedie_p, sigma_O, sigma_E, sigma_Z, rho)
}

print_header <- function(x) {
  info <- print_model_info(x)
  cat(info$title)
  cat(info$formula)
  cat(info$mesh)
  cat(info$time)
  cat(info$data)
  cat(info$overall_family)
}

print_one_model <- function(x, m = 1) {
  info <- print_model_info(x)
  main <- print_main_effects(x, m = m)
  smooth <- print_smooth_effects(x, m = m)
  iid_re <- print_iid_re(x, m = m)
  tv <- print_time_varying(x, m = m)
  range <- print_range(x, m = m)
  other <- print_other_parameters(x, m = m)

  if (m == 1) cat(info$family1, "\n")
  if (m == 2) cat(info$family2, "\n")

  print(rbind(main, smooth$smooth_effects))
  cat("\n")

  if (!is.null(smooth$smooth_sds)) {
    cat("Smooth terms:\n")
    print(smooth$smooth_sds)
    cat("\n")
  }

  if (!is.null(iid_re)) {
    cat("Random intercepts:\n")
    print(iid_re)
    cat("\n")
  }

  if (!is.null(tv)) {
    cat("Time-varying parameters:\n")
    print(tv)
    cat("\n")
  }

  cat(other$phi)
  cat(other$tweedie_p)
  cat(other$rho)
  cat(range)
  cat(other$sigma_O)
  cat(other$sigma_Z)
  cat(other$sigma_E)
}
print_footer <- function(x) {
  info <- print_model_info(x)
  cat(info$criterion)

  cat("\nSee ?tidy.sdmTMB to extract these values as a data frame.\n")

  # suppressWarnings(suppressMessages(ok <- sanity(x)$all_ok))
  # if (!isTRUE(ok)) {
    # cat("Possible issues detected; printing output of sanity():\n\n")
    # sanity(x)
  # }
}

#' @export
#' @import methods
print.sdmTMB <- function(x, ...) {

  # or x$tmb_obj$retape()!?
  sink(tempfile())
  tmp <- x$tmb_obj$fn(x$tmb_obj$par) # FIXME needed?
  lp <- x$tmb_obj$env$last.par.best
  r <- x$tmb_obj$report(lp)
  sink()

  delta <- isTRUE(x$family$delta)
  print_header(x)
  if (delta) cat("\nDelta/hurdle model 1: -----------------------------------\n")
  print_one_model(x, 1)
  if (delta) {
    cat("\nDelta/hurdle model 2: -----------------------------------\n")
    print_one_model(x, 2)
  }
  if (delta) cat("\n")
  print_footer(x)
}

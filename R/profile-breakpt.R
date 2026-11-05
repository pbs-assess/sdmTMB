#' Estimate a break-point cutpoint by profile likelihood
#'
#' Fits a model with a `breakpt()` term (see the threshold models section of
#' [sdmTMB()]) over a grid of fixed cutpoints, finds the cutpoint with the
#' lowest negative log likelihood (NLL), and optionally refits the model with
#' the cutpoint fixed at that value. Use this in place of estimating the
#' cutpoint directly with [sdmTMB()].
#'
#' @details
#' The NLL of a `breakpt()` model has a kink in the cutpoint at every
#' observation, so an optimizer that estimates the cutpoint with the other
#' parameters can stop at the wrong cutpoint while reporting convergence, and
#' its standard error for the cutpoint isn't reliable. With the cutpoint held
#' fixed, the likelihood is smooth in the other parameters.
#'
#' The profile works in three steps:
#'
#' 1. Fit the model at each cutpoint in `cuts`, working outward from the median
#'    of the covariate.
#' 2. Wherever the profile is within `within` NLL units of its best, add
#'    cutpoints until neighbouring cutpoints are at most `tol_coarse` apart.
#'    Next to the best cutpoint, and on either side of each point where the
#'    profile crosses the interval cutoff, continue until they are at most
#'    `tol` apart.
#' 3. Search for the minimum between the best cutpoint's two neighbours with
#'    [stats::optimize()].
#'
#' Each profile fit starts from the nearest finished fit and skips the extra
#' optimization steps and [TMB::sdreport()] for speed.
#' Only the final fit runs those.
#'
#' The confidence interval is every cutpoint within `qchisq(level, 1) / 2` NLL units of
#' the best, with the ends interpolated between grid points. If the profile
#' is still inside the cutoff at the end of the grid, that end is `NA`.
#'
#' @param formula,... As in [sdmTMB()]. `formula` must include one `breakpt()`
#'   term. Pass `data`, `mesh`, `family`, and other [sdmTMB()] arguments in
#'   `...`.
#' @param control Optimization control options from [sdmTMBcontrol()]. Used as
#'   given for the final fit, except that `multiphase` is turned off, at least
#'   one Newton loop runs, and `start` and `map` are set to fix the cutpoint.
#'   Don't map `b_threshold` yourself.
#' @param cuts Cutpoints for the first pass, on the scale of the covariate. By
#'   default, quantiles of the covariate from 1% to 99%, with extra points in
#'   the tails, plus 8 evenly spaced points over the same range. The outer 1%
#'   at each end is left out because so few observations fall beyond such a
#'   cutpoint.
#' @param within Refine wherever the profile is within this many NLL units of
#'   its best.
#' @param tol Final spacing of cutpoints next to the best cutpoint and around
#'   the ends of the interval. Defaults to 0.025 times the standard deviation
#'   of the covariate.
#' @param tol_coarse Spacing of cutpoints elsewhere within `within`, fine
#'   enough to find a second minimum. Defaults to 4 times `tol`.
#' @param max_fits Stop refining after this many profile fits. The search in
#'   step 3 runs regardless.
#' @param level Confidence level for the profile interval.
#' @param final Fit the final model at the best cutpoint? If `FALSE`, returns
#'   the profile only.
#'
#' @return A list of class `breakpt_profile` with:
#' * `profile`: a data frame with one row per profile fit: the cutpoint
#'   (`cut`), its NLL (`nll`, `NA` if the fit failed), NLL above the best
#'   (`dnll`), the slope, the maximum absolute gradient, the convergence code,
#'   and the share of observations below the cutpoint.
#' * `estimate`: a one-row data frame with the best cutpoint (`cut`), slope,
#'   NLL, interval (`lwr`, `upr`), AIC, numbers of observations below and
#'   above the cutpoint, and flags for checking the result (see Warnings).
#' * `fit`: the [sdmTMB()] fit with the cutpoint fixed at `estimate$cut`, or
#'   `NULL` if `final = FALSE`.
#'
#' Use `estimate$aic`, not `AIC(estimate$fit)`: `fit` doesn't count the fixed cutpoint
#' as a parameter. Standard errors from `fit` don't include uncertainty in the
#' cutpoint.
#'
#' @section Warnings:
#' The function warns, and sets `estimate$open_interval` or
#' `estimate$few_below`, when the interval reaches the end of the grid or when
#' fewer than 15 observations (or 5%, if that's more) fall below the best
#' cutpoint. Either means the data say little about the cutpoint. Also check
#' `estimate$pieces` (the interval isn't continuous if greater than 1),
#' `estimate$near_edge`, `estimate$refine_complete`, and
#' `estimate$max_grad_near` (the largest gradient among fits within `within`
#' NLL units of the best).
#'
#' Delta models and models with more than one threshold term aren't
#' supported.
#'
#' @seealso [plot.breakpt_profile()]
#' @export
#' @examples
#' # Simulate a hockey-stick relationship with a cutpoint at 0.4:
#' set.seed(1)
#' d <- data.frame(x = runif(200, -2, 2))
#' d$y <- 1 + 1.5 * pmin(d$x, 0.4) + rnorm(200, sd = 0.5)
#'
#' p <- profile_breakpt(y ~ breakpt(x), data = d, spatial = "off")
#' p
#' if (ggplot2_installed()) plot(p)
#'
#' p$estimate[, c("cut", "lwr", "upr", "aic")]
#' head(p$profile)
#'
#' # The final fit, with the cutpoint fixed at its best value:
#' p$fit
#'
#' # The fitted hockey stick on the data, with the best cutpoint (dashed) and
#' # its 95% interval (dotted):
#' nd <- data.frame(x = seq(min(d$x), max(d$x), length.out = 200))
#' nd$est <- predict(p$fit, newdata = nd)$est
#' plot(d$x, d$y, col = "grey60", xlab = "x", ylab = "y")
#' lines(nd$x, nd$est, lwd = 2)
#' abline(v = p$estimate$cut, lty = 2)
#' abline(v = c(p$estimate$lwr, p$estimate$upr), lty = 3)
profile_breakpt <- function(formula, ..., control = sdmTMBcontrol(), cuts = NULL,
                            within = 4, tol = NULL, tol_coarse = NULL,
                            max_fits = 100L, level = 0.95, final = TRUE) {
  if (!is.null(control$map$b_threshold)) {
    cli_abort("Don't map `b_threshold` in `control`; the profile sets it.")
  }
  control_build <- control
  control_build$normalize <- FALSE # normalizing only saves time inside the inner problem
  built <- suppressMessages(sdmTMB(formula, ..., control = control_build, do_fit = FALSE))
  obj <- built$tmb_obj
  if (!identical(as.integer(built$tmb_data$threshold_func), 1L)) {
    cli_abort("`formula` needs a `breakpt()` term.")
  }
  ib <- which(names(obj$par) == "b_threshold")
  if (length(ib) != 2L) cli_abort("Delta models aren't supported.")
  i_slope <- ib[1]
  i_cut <- ib[2]
  x <- as.numeric(built$tmb_data$X_threshold)
  if (is.null(tol)) tol <- 0.025 * stats::sd(x)
  if (is.null(tol_coarse)) tol_coarse <- 4 * tol
  # otherwise refinement may never stop
  if (!(tol > 0 && tol_coarse > 0 && is.finite(max_fits))) {
    cli_abort("`tol` and `tol_coarse` must be positive and `max_fits` finite.")
  }
  crit <- stats::qchisq(level, 1) / 2
  if (is.null(cuts)) {
    probs <- c(0.01, 0.015, 0.025, seq(0.05, 0.95, length.out = 15), 0.975, 0.985, 0.99)
    ends <- stats::quantile(x, c(0.01, 0.99), names = FALSE)
    cuts <- c(stats::quantile(x, probs, names = FALSE), seq(ends[1], ends[2], length.out = 8))
  }
  cuts <- sort(unique(cuts))
  if (length(cuts) < 3L) cli_abort("Give at least 3 `cuts`.")

  fit_at <- function(ct, start_par, start_last) {
    start_par[i_cut] <- ct
    free <- setdiff(seq_along(start_par), i_cut)
    full <- start_par
    obj$env$last.par.best <- start_last
    obj$env$value.best <- Inf
    fn <- function(p) {
      full[free] <- p
      v <- obj$fn(full)
      if (is.finite(v)) v else 1e10
    }
    gr <- function(p) {
      full[free] <- p
      g <- obj$gr(full)[free]
      g[!is.finite(g)] <- 0
      g
    }
    o <- tryCatch(
      suppressWarnings(stats::nlminb(start_par[free], fn, gr,
        control = list(eval.max = 2000L, iter.max = 1000L))),
      error = function(e) NULL
    )
    # 1e10 means nlminb() never reached a finite NLL
    if (is.null(o) || !is.finite(o$objective) || o$objective >= 1e10) return(NULL)
    full[free] <- o$par
    list(cut = ct, nll = o$objective, par = full, last = obj$env$last.par.best,
      max_grad = max(abs(gr(o$par))), convergence = o$convergence)
  }

  done <- list()
  add_fit <- function(ct) {
    if (length(done)) {
      have <- vapply(done, `[[`, 0, "cut")
      n <- done[[which.min(abs(have - ct))]]
      r <- fit_at(ct, n$par, n$last)
    } else {
      r <- fit_at(ct, obj$par, obj$env$last.par.best)
    }
    if (is.null(r)) {
      r <- list(cut = ct, nll = NA_real_, par = NULL, last = NULL,
        max_grad = NA_real_, convergence = NA_integer_)
    }
    done[[length(done) + 1L]] <<- r
    r
  }
  # a failed fit stays in `done` so it isn't retried, but can't seed a neighbour
  add_fit_safe <- function(ct) {
    keep <- done
    done <<- Filter(function(z) !is.null(z$par), done)
    r <- add_fit(ct)
    done <<- c(keep, list(r))
    r
  }

  # First pass: out from the median in both directions, so each fit starts from its neighbour
  first <- which.min(abs(cuts - stats::median(x)))
  add_fit_safe(cuts[first])
  if (is.null(done[[1]]$par)) cli_abort("The first profile fit failed. Try other `cuts` or check the model.")
  for (ct in cuts[seq_along(cuts) > first]) add_fit_safe(ct)
  for (ct in rev(cuts[seq_along(cuts) < first])) add_fit_safe(ct)

  table_now <- function() {
    d <- data.frame(cut = vapply(done, `[[`, 0, "cut"), nll = vapply(done, `[[`, 0, "nll"))
    d[order(d$cut), ]
  }
  # Refine: bisect gaps near the best until they are at most `tol_coarse`, and the gaps that
  # set the answer (next to the best, and across each crossing of the interval cutoff) until
  # they are at most `tol`. Every pass adds a fit, so this runs at most `max_fits` times.
  repeat {
    d <- table_now()
    dn <- d$nll - min(d$nll, na.rm = TRUE)
    near <- !is.na(dn) & dn < within
    gap <- diff(d$cut)
    left <- dn[-nrow(d)]
    right <- dn[-1]
    i_best <- which.min(dn)
    crosses <- !is.na(left) & !is.na(right) & (left < crit) != (right < crit)
    key <- crosses | seq_along(gap) %in% c(i_best - 1L, i_best)
    todo <- which(((near[-1] | near[-nrow(d)]) & gap > tol_coarse) | (key & gap > tol))
    if (!length(todo) || length(done) >= max_fits) break
    new <- (d$cut[todo] + d$cut[todo + 1L]) / 2
    # closest to the current best first, so warm starts move outward from it
    new <- new[order(abs(new - d$cut[i_best]))]
    for (ct in new) if (length(done) < max_fits) add_fit_safe(ct)
  }
  refine_complete <- !length(todo)

  # Search between the best cutpoint's neighbours. The profile is continuous, even where it
  # has kinks, so a 1-D search is safe. Every evaluation is kept as a profile fit.
  d <- table_now()
  d <- d[!is.na(d$nll), ]
  i_best <- which.min(d$nll)
  bracket <- d$cut[c(max(1L, i_best - 1L), min(nrow(d), i_best + 1L))]
  if (diff(bracket) > 0) {
    stats::optimize(function(ct) {
      have <- vapply(done, `[[`, 0, "cut")
      r <- if (ct %in% have) done[[match(ct, have)]] else add_fit_safe(ct) # optimize() can repeat a point
      if (is.na(r$nll)) 1e10 else r$nll
    }, bracket, tol = tol / 20)
  }

  done <- done[order(vapply(done, `[[`, 0, "cut"))]
  p <- data.frame(
    cut = vapply(done, `[[`, 0, "cut"),
    nll = vapply(done, `[[`, 0, "nll"),
    slope = vapply(done, function(z) if (is.null(z$par)) NA_real_ else z$par[[i_slope]], 0),
    max_grad = vapply(done, `[[`, 0, "max_grad"),
    convergence = vapply(done, function(z) as.integer(z$convergence), 0L)
  )
  p$dnll <- p$nll - min(p$nll, na.rm = TRUE)
  p$prop_below <- vapply(p$cut, function(ct) mean(x < ct), 0)
  best_i <- which.min(p$nll)
  best <- done[[best_i]]

  # Interval: outermost crossings of the cutoff, interpolated between grid points
  q <- p[!is.na(p$nll), ]
  inside <- q$dnll < crit
  crossing <- function(i, j) {
    q$cut[i] + (crit - q$dnll[i]) / (q$dnll[j] - q$dnll[i]) * (q$cut[j] - q$cut[i])
  }
  il <- min(which(inside))
  iu <- max(which(inside))
  lo <- if (il > 1L) crossing(il, il - 1L) else NA_real_ # NA: still inside at the end of the grid
  hi <- if (iu < nrow(q)) crossing(iu, iu + 1L) else NA_real_

  k <- length(obj$par) # includes the cutpoint
  n_below <- sum(x < best$cut)
  estimate <- data.frame(
    cut = best$cut, slope = p$slope[best_i], nll = best$nll, lwr = lo, upr = hi,
    pieces = sum(diff(c(FALSE, inside)) == 1L),
    k = k, aic = 2 * best$nll + 2 * k, prop_below = p$prop_below[best_i],
    n_below = n_below, n_above = length(x) - n_below,
    near_edge = best_i <= 2L || best_i >= nrow(p) - 1L, refine_complete = refine_complete,
    n_fits = nrow(p), n_failed = sum(is.na(p$nll)),
    # large gradients matter only near the best; far away they can't change the result
    max_grad_near = max(p$max_grad[!is.na(p$dnll) & p$dnll < within], na.rm = TRUE),
    open_interval = is.na(lo) || is.na(hi),
    few_below = n_below < max(15, 0.05 * length(x))
  )
  if (estimate$open_interval) {
    cli_warn(paste0("The ", 100 * level, "% interval for the cutpoint reaches the end of ",
      "the grid, so the data don't rule out a cutpoint at the edge of their range."))
  }
  if (estimate$few_below) {
    cli_warn(paste0("Only ", n_below, " observations fall below the best cutpoint, so ",
      "the cutpoint and slope may be poorly identified."))
  }

  fit <- NULL
  if (final) {
    ctl <- control
    ctl$multiphase <- FALSE
    ctl$newton_loops <- max(1L, ctl$newton_loops)
    ctl$start <- utils::modifyList(as.list(ctl$start), obj$env$parList(best$par, best$last))
    ctl$map <- c(ctl$map, list(b_threshold = factor(c(1L, NA))))
    fit <- suppressMessages(sdmTMB(formula, ..., control = ctl))
    estimate$final_nll <- fit$model$objective
    estimate$nll <- min(estimate$nll, fit$model$objective)
    estimate$aic <- 2 * estimate$nll + 2 * k
  }
  structure(list(profile = p, estimate = estimate, fit = fit, level = level,
    within = within, tol = tol, tol_coarse = tol_coarse), class = "breakpt_profile")
}

#' @export
print.breakpt_profile <- function(x, ...) {
  e <- x$estimate
  cat(sprintf("Break-point profile: %d fits (%d failed)\n", e$n_fits, e$n_failed))
  cat(sprintf("Cutpoint %.3f, %g%% interval %.3f to %.3f\n",
    e$cut, 100 * x$level, e$lwr, e$upr))
  cat(sprintf("Slope %.3f, NLL %.3f, AIC %.2f (k = %d, counting the cutpoint)\n",
    e$slope, e$nll, e$aic, e$k))
  cat(sprintf("Observations below the cutpoint: %d (%.1f%%)\n", e$n_below, 100 * e$prop_below))
  if (e$open_interval) cat("Note: the interval reaches the end of the grid.\n")
  if (e$few_below) cat("Note: few observations fall below the cutpoint.\n")
  if (e$pieces > 1L) cat("Note: the interval is not continuous; see `$profile`.\n")
  if (e$near_edge) cat("Note: the best cutpoint is at the edge of the grid.\n")
  if (e$max_grad_near > 0.1) cat("Note: a profile fit near the best has a large gradient; see `$profile`.\n")
  if (!e$refine_complete) cat("Note: refinement stopped at `max_fits`.\n")
  invisible(x)
}

#' Plot a break-point profile
#'
#' Plots the profile likelihood from [profile_breakpt()]: the NLL above its
#' minimum at each cutpoint tried, the cutoff for the confidence interval, the
#' interval itself, and the best cutpoint. Needs the ggplot2 package.
#'
#' @param x Output from [profile_breakpt()].
#' @param ylim Range of the vertical axis. The profile is cut off above it, so
#'   the region near the minimum stays readable.
#' @param ... Not used.
#'
#' @return A ggplot object. The shaded band is the confidence interval; if the
#'   interval is open at an end, the band runs to the end of the grid. Tick
#'   marks along the bottom show the cutpoints tried.
#' @export
#' @examplesIf ggplot2_installed()
#' set.seed(1)
#' d <- data.frame(x = runif(200, -2, 2))
#' d$y <- 1 + 1.5 * pmin(d$x, 0.4) + rnorm(200, sd = 0.5)
#' p <- profile_breakpt(y ~ breakpt(x), data = d, spatial = "off", final = FALSE)
#' plot(p)
#' plot(p, ylim = c(0, 25))
plot.breakpt_profile <- function(x, ylim = c(0, 10), ...) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    cli_abort("ggplot2 must be installed to use this function.")
  }
  col <- "#1b9e6e"
  ink <- "grey25"
  p <- x$profile[!is.na(x$profile$dnll), ]
  e <- x$estimate
  crit <- stats::qchisq(x$level, 1) / 2
  band <- data.frame(
    lwr = if (is.na(e$lwr)) min(p$cut) else e$lwr,
    upr = if (is.na(e$upr)) max(p$cut) else e$upr
  )
  best <- data.frame(cut = e$cut, dnll = 0)
  # the cutoff label goes at the end of the grid farther from the best cutpoint, clear of the curve
  right <- e$cut < mean(range(p$cut))
  ggplot2::ggplot(p, ggplot2::aes(.data$cut, .data$dnll)) +
    ggplot2::geom_rect(data = band, ggplot2::aes(xmin = .data$lwr, xmax = .data$upr,
      ymin = -Inf, ymax = Inf), inherit.aes = FALSE, fill = col, alpha = 0.12) +
    ggplot2::geom_hline(yintercept = crit, linetype = "22", colour = ink, linewidth = 0.4) +
    ggplot2::annotate("text", x = if (right) max(p$cut) else min(p$cut), y = crit,
      label = paste0(100 * x$level, "% cutoff"), hjust = if (right) 1 else 0, vjust = -0.6,
      size = 3.2, colour = ink) +
    # the profile is below the cutoff at the best cutpoint, so its label above the cutoff is clear
    ggplot2::annotate("text", x = e$cut, y = crit, label = sprintf("%.2f", e$cut),
      vjust = -0.6, size = 3.2, colour = ink) +
    ggplot2::geom_line(colour = col, linewidth = 0.9) +
    ggplot2::geom_rug(sides = "b", colour = ink, linewidth = 0.3,
      length = ggplot2::unit(0.03, "npc")) +
    ggplot2::geom_point(data = best, size = 3.2, shape = 21, fill = "white", colour = col,
      stroke = 1.3) +
    ggplot2::coord_cartesian(ylim = ylim) +
    ggplot2::labs(x = "Cutpoint", y = "NLL above the minimum") +
    ggplot2::theme_minimal() +
    ggplot2::theme(panel.grid.minor = ggplot2::element_blank())
}

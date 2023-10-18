#' Run extra optimization on an already fitted object
#'
#' @description
#' `r lifecycle::badge("experimental")`
#'
#' @param object An object from [sdmTMB()].
#' @param nlminb_loops How many extra times to run [stats::nlminb()]
#'   optimization. Sometimes restarting the optimizer at the previous best
#'   values aids convergence.
#' @param newton_loops How many extra Newton optimization loops to try with
#'   [stats::optimHess()]. Sometimes aids convergence.
#'
#' @return An updated model fit of class `sdmTMB`.
#' @export
#'
#' @examples
#' # Run extra optimization steps to help convergence:
#' # (Not typically needed)
#' fit <- sdmTMB(density ~ 0 + poly(depth, 2) + as.factor(year),
#'   data = pcod_2011, mesh = pcod_mesh_2011, family = tweedie())
#' fit_1 <- run_extra_optimization(fit, newton_loops = 1)
#' max(fit$gradients)
#' max(fit_1$gradients)
run_extra_optimization <- function(object,
  nlminb_loops = 0,
  newton_loops = 1) {
  # FIXME: DRY; use this function in sdmTMB()

  new_obj <- object
  old_par <- object$model$par
  new_obj$tmb_obj$fn(old_par) # initialize the TMB object
  tmb_opt <- new_obj$model
  for (i in seq_len(nlminb_loops)) {
    temp <- new_obj$model[c("iterations", "evaluations")]
    tmb_opt <- stats::nlminb(
      start = tmb_opt$par,
      objective = new_obj$tmb_obj$fn,
      gradient = new_obj$tmb_obj$gr,
      control = object$nlminb_control,
      lower = object$lower,
      upper = object$upper
    )
    tmb_opt[["iterations"]] <- tmb_opt[["iterations"]] + temp[["iterations"]]
    tmb_opt[["evaluations"]] <- tmb_opt[["evaluations"]] + temp[["evaluations"]]
  }
  for (i in seq_len(newton_loops)) {
    g <- as.numeric(new_obj$tmb_obj$gr(tmb_opt$par))
    h <- stats::optimHess(
      tmb_opt$par,
      fn = new_obj$tmb_obj$fn,
      gr = new_obj$tmb_obj$gr,
      lower = object$lower,
      upper = object$upper
    )
    tmb_opt$par <- tmb_opt$par - solve(h, g)
    tmb_opt$objective <- new_obj$tmb_obj$fn(tmb_opt$par)
  }
  new_obj$model <- tmb_opt
  new_obj$sd_report <- TMB::sdreport(new_obj$tmb_obj,
    getJointPrecision = "jointPrecision" %in% names(object$sd_report))
  conv <- get_convergence_diagnostics(new_obj$sd_report)
  new_obj$gradients <- conv$final_grads
  new_obj$bad_eig <- conv$bad_eig
  new_obj
}

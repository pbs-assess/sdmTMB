#' Run extra optimization on an already fitted object
#'
#' @param object An object from [sdmTMB()].
#' @param nlminb_loops How many extra times to run [stats::nlminb()]
#'   optimization. Sometimes restarting the optimizer at the previous best
#'   values aids convergence.
#' @param newton_steps How many extra Newton optimization steps to try with
#'   [stats::optimHess()]. Sometimes aids convergence.
#' @param control Optimization control options. See [sdmTMBcontrol()].
#'
#' @return An updated model fit of class `sdmTMB`.
#' @export
#'
#' @examples
#' # See ?sdmTMB
run_extra_optimization <- function(object,
  nlminb_loops = 1,
  newton_steps = 1,
  control = sdmTMBcontrol()) {
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
      control = control
    )
    tmb_opt[["iterations"]] <- tmb_opt[["iterations"]] + temp[["iterations"]]
    tmb_opt[["evaluations"]] <- tmb_opt[["evaluations"]] + temp[["evaluations"]]
  }
  for (i in seq_len(newton_steps)) {
    g <- as.numeric(new_obj$tmb_obj$gr(tmb_opt$par))
    h <- stats::optimHess(
      tmb_opt$par,
      fn = new_obj$tmb_obj$fn,
      gr = new_obj$tmb_obj$gr
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

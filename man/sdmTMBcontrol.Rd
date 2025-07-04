% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{sdmTMBcontrol}
\alias{sdmTMBcontrol}
\title{Optimization control options}
\usage{
sdmTMBcontrol(
  eval.max = 2000L,
  iter.max = 1000L,
  normalize = FALSE,
  nlminb_loops = 1L,
  newton_loops = 1L,
  mgcv = deprecated(),
  quadratic_roots = FALSE,
  start = NULL,
  map_rf = deprecated(),
  map = NULL,
  lower = NULL,
  upper = NULL,
  censored_upper = NULL,
  multiphase = TRUE,
  profile = FALSE,
  get_joint_precision = TRUE,
  parallel = getOption("sdmTMB.cores", 1L),
  suppress_nlminb_warnings = TRUE,
  ...
)
}
\arguments{
\item{eval.max}{Maximum number of evaluations of the objective function
allowed.}

\item{iter.max}{Maximum number of iterations allowed.}

\item{normalize}{Logical: use \code{\link[TMB:normalize]{TMB::normalize()}} to normalize the process
likelihood using the Laplace approximation? Can result in a substantial
speed boost in some cases. This used to default to \code{FALSE} prior to
May 2021. Currently not working for models fit with REML or random intercepts.}

\item{nlminb_loops}{How many times to run \code{\link[stats:nlminb]{stats::nlminb()}} optimization.
Sometimes restarting the optimizer at the previous best values aids
convergence. If the maximum gradient is still too large,
try increasing this to \code{2}.}

\item{newton_loops}{How many Newton optimization steps to try after running
\code{\link[stats:nlminb]{stats::nlminb()}}. This sometimes aids convergence by further reducing the
log-likelihood gradient with respect to the fixed effects. This calculates
the Hessian at the current MLE with \code{\link[stats:optim]{stats::optimHess()}} using a
finite-difference approach and uses this to update the fixed effect
estimates.}

\item{mgcv}{\strong{Deprecated} Parse the formula with \code{\link[mgcv:gam]{mgcv::gam()}}?}

\item{quadratic_roots}{Experimental feature for internal use right now; may
be moved to a branch. Logical: should quadratic roots be calculated? Note:
on the sdmTMB side, the first two coefficients are used to generate the
quadratic parameters. This means that if you want to generate a quadratic
profile for depth, and depth and depth^2 are part of your formula, you need
to make sure these are listed first and that an intercept isn't included.
For example, \code{formula = cpue ~ 0 + depth + depth2 + as.factor(year)}.}

\item{start}{A named list specifying the starting values for parameters. You
can see the necessary structure by fitting the model once and inspecting
\code{your_model$tmb_obj$env$parList()}. Elements of \code{start} that are specified
will replace the default starting values.}

\item{map_rf}{\strong{Deprecated} use \verb{spatial = 'off', spatiotemporal = 'off'} in
\code{\link[=sdmTMB]{sdmTMB()}}.}

\item{map}{A named list with factor \code{NA}s specifying parameter values that
should be fixed at a constant value. See the documentation in
\code{\link[TMB:MakeADFun]{TMB::MakeADFun()}}. This should usually be used with \code{start} to specify the
fixed value.}

\item{lower}{An optional named list of lower bounds within the optimization.
Parameter vectors with the same name (e.g., \code{b_j} or \code{ln_kappa} in some
cases) can be specified as a numeric vector. E.g.
\code{lower = list(b_j = c(-5, -5))}. Note that \code{\link[stats:optim]{stats::optimHess()}} does not
implement lower and upper bounds, so you must set \code{newton_loops = 0} if
setting limits.}

\item{upper}{An optional named list of upper bounds within the optimization.}

\item{censored_upper}{An optional vector of upper bounds for
\code{\link[=sdmTMBcontrol]{sdmTMBcontrol()}}. Values of \code{NA} indicate an unbounded right-censored to
distribution, values greater that the observation indicate and upper bound,
and values equal to the observation indicate no censoring.}

\item{multiphase}{Logical: estimate the fixed and random effects in phases?
Phases are usually faster and more stable.}

\item{profile}{Logical: should population-level/fixed effects be profiled
out of the likelihood? These are then appended to the random effects
vector without the Laplace approximation. See \code{\link[TMB:MakeADFun]{TMB::MakeADFun()}}. \emph{This
can dramatically speed up model fit if there are many fixed effects but is
experimental at this stage.}}

\item{get_joint_precision}{Logical. Passed to \code{getJointPrecision} in
\code{\link[TMB:sdreport]{TMB::sdreport()}}. Must be \code{TRUE} to use simulation-based methods in
\code{\link[=predict.sdmTMB]{predict.sdmTMB()}} or \verb{[get_index_sims()]}. If not needed, setting this
\code{FALSE} will reduce object size.}

\item{parallel}{Argument currently ignored. For parallel processing with 3
cores, as an example, use \code{TMB::openmp(n = 3, DLL = "sdmTMB")}. But be
careful, because it's not always faster with more cores and there is
definitely an upper limit.}

\item{suppress_nlminb_warnings}{Suppress uninformative warnings
from \code{\link[stats:nlminb]{stats::nlminb()}} arising when a function evaluation is \code{NA}, which
are then replaced with \code{Inf} and avoided during estimation?}

\item{...}{Anything else. See the 'Control parameters' section of
\code{\link[stats:nlminb]{stats::nlminb()}}.}
}
\value{
A list of control arguments
}
\description{
\code{\link[=sdmTMB]{sdmTMB()}} and \code{\link[stats:nlminb]{stats::nlminb()}} control options.
}
\details{
Usually used within \code{\link[=sdmTMB]{sdmTMB()}}. For example:

\if{html}{\out{<div class="sourceCode">}}\preformatted{sdmTMB(..., control = sdmTMBcontrol(newton_loops = 2))
}\if{html}{\out{</div>}}
}
\examples{
sdmTMBcontrol()
}

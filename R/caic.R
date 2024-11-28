
#' @title Calculate conditional AIC
#'
#' @description
#' Calculates the conditional Akaike Information criterion (cAIC).
#'
#' @param object Output from \code{\link{sdmTMB}}
#' @param what Whether to return the cAIC or the effective degrees of freedom
#'        (EDF) for each group of random effects.
#'
#' @details
#' cAIC is designed to optimize the expected out-of-sample predictive
#' performance for new data that share the same random effects as the
#' in-sample (fitted) data, e.g., spatial interpolation.  In this sense,
#' it should be a fast approximation to optimizing the model structure
#' based on k-fold crossvalidation.
#' By contrast, \code{AIC} calculates the
#' marginal Akaike Information Criterion, which is designed to optimize
#' expected predictive performance for new data that have new random effects,
#' e.g., extrapolation, or inference about generative parameters.
#'
#' cAIC also calculates as a byproduct the effective degrees of freedom,
#' i.e., the number of fixed effects that would have an equivalent impact on
#' model flexibility as a given random effect.
#'
#' Both cAIC and EDF are calculated using Eq. 6 of Zheng Cadigan Thorson 2024.
#'
#' Note that, for models that include profiled fixed effects, these profiles
#' are turned off.
#'
#' @return
#' Either the cAIC, or the effective degrees of freedom (EDF) by group
#' of random effects
#'
#' @references
#'
#' **Deriving the general approximation to cAIC used here**
#'
#' Zheng, N., Cadigan, N., & Thorson, J. T. (2024).
#' A note on numerical evaluation of conditional Akaike information for
#' nonlinear mixed-effects models (arXiv:2411.14185). arXiv.
#' \doi{10.48550/arXiv.2411.14185}
#'
#' **The utility of EDF to diagnose hierarchical model behavior**
#'
#' Thorson, J. T. (2024). Measuring complexity for hierarchical
#' models using effective degrees of freedom. Ecology,
#' 105(7), e4327 \doi{10.1002/ecy.4327}
#'
#' @export
CAIC.sdmTMB <-
function( object,
          what = c("CAIC","EDF") ){

  what = match.arg(what)
  require(Matrix)
  tmb_data = object$tmb_data

  # Make sure profile = NULL
  if( is.null(object$control$profile) ){
    obj = object$tmb_obj
  }else{
    obj = TMB::MakeADFun( data = tmb_data,
                         parameters = object$parlist,
                         map = object$tmb_map,
                         random = object$tmb_random,
                         DLL = "sdmTMB",
                         profile = NULL )
  }

  # Make obj_new
  tmb_data$weights_i[] = 0
  obj_new = TMB::MakeADFun( data = tmb_data,
                      parameters = object$parlist,
                      map = object$tmb_map,
                      random = object$tmb_random,
                      DLL = "sdmTMB",
                      profile = NULL )

  #
  par = obj$env$parList()
  parDataMode <- obj$env$last.par
  indx = obj$env$lrandom()
  q = length(indx)
  p = length(object$model$par)

  ## use - for Hess because model returns negative loglikelihood;
  #cov_Psi_inv = -Hess_new[indx,indx]; ## this is the marginal prec mat of REs;
  Hess_new = -Matrix(obj_new$env$f(parDataMode,order=1,type="ADGrad"),sparse = TRUE)
  Hess_new = Hess_new[indx,indx]

  ## Joint hessian etc
  Hess = -Matrix(obj$env$f(parDataMode,order=1,type="ADGrad"),sparse = TRUE)
  Hess = Hess[indx,indx]
  negEDF = diag(solve(Hess, Hess_new))

  if(what=="CAIC"){
    jnll = obj$env$f(parDataMode)
    cnll = jnll - obj_new$env$f(parDataMode)
    cAIC = 2*cnll + 2*(p+q) - 2*sum(negEDF)
    return(cAIC)
  }
  if(what=="EDF"){
    # Figure out group for each random-effect coefficient
    group = factor(names(object$last.par.best[obj$env$random]))
    # Calculate total EDF by group
    EDF = tapply(negEDF,INDEX=group,FUN=length) - tapply(negEDF,INDEX=group,FUN=sum)
    return(EDF)
  }
}

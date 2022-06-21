// PC Prior only
#include <TMB.hpp>
template<class Type>
  Type objective_function<Type>::operator() ()
{
  // these are data inputs / hyperparameters for PC prior
  DATA_SCALAR(matern_range);
  DATA_SCALAR(range_prob);
  DATA_SCALAR(matern_SD);
  DATA_SCALAR(SD_prob);
  DATA_INTEGER(include_jacob);

  // parameters
  PARAMETER(ln_tau_O);    // spatial process
  PARAMETER(ln_kappa);    // Matern parameter

  Type range = sqrt(Type(8.)) / exp(ln_kappa); // range as derived parameter

  // This is taken from utils.h in sdmTMB
  Type d = 2.;  // dimension
  Type dhalf = d / 2.;
  Type lam1 = -log(range_prob) * pow(matern_range, dhalf);
  Type lam2 = -log(SD_prob) / matern_SD;
  Type sigma = 1 / sqrt(Type(4.) * M_PI * exp(Type(2.) * ln_tau_O + Type(2.) * ln_kappa));

  Type range_ll = log(dhalf) + log(lam1) + log(pow(range, -1. - dhalf)) -
    lam1 * pow(range, -dhalf);
  Type sigma_ll = log(lam2) - lam2 * sigma;
  Type penalty = sigma_ll + range_ll;

  Type(jnll) = 0.;
  // add prior penalty
  jnll -= penalty;

  // add Jacobian adjustments?
  if (include_jacob) {
    jnll -= log(sqrt(8.0)) - log(pow(range,2.0)); // P(sigma)
    Type C = sqrt(exp(lgamma(1.0 + dhalf)) * pow(4*M_PI, dhalf));
    jnll -= log(C) + ln_kappa;
  }

  // from calc_rf
  Type sigma_O = 1 / sqrt(Type(4.) * M_PI * exp(Type(2.) * ln_tau_O + Type(2.) * ln_kappa));
  REPORT(range);
  ADREPORT(range);
  REPORT(sigma_O);
  ADREPORT(sigma_O);
  return jnll;
}

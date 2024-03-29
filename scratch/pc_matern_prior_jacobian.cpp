// Simple spatial model.
#include <TMB.hpp>
template<class Type>
  Type objective_function<Type>::operator() ()
{
  using namespace R_inla;
  using namespace density;
  using namespace Eigen;

  // these are data inputs / hyperparameters for PC prior
  DATA_SCALAR(matern_range);
  DATA_SCALAR(range_prob);
  DATA_SCALAR(matern_SD);
  DATA_SCALAR(SD_prob);
  DATA_STRUCT(spde, spde_t);
  DATA_SPARSE_MATRIX(A_st); // INLA 'A' projection matrix for unique stations
  DATA_IVECTOR(A_spatial_index); // Vector of stations to match up A_st output
  DATA_VECTOR(y);

  // parameters
  PARAMETER(B0); // intercept
  PARAMETER(ln_tau_O);    // spatial process
  PARAMETER(ln_kappa);    // Matern parameter
  PARAMETER(ln_sigma);    // observation sd
  PARAMETER_ARRAY(omega_s); // random spatial component

  int n_i = y.rows();   // number of observations

  Type range = sqrt(Type(8.)) / exp(ln_kappa); // range as derived parameter
  Eigen::SparseMatrix<Type> Q_s; // Precision matrix
  Q_s = R_inla::Q_spde(spde, exp(ln_kappa));
  Type jnll = 0;
  jnll += SCALE(GMRF(Q_s, true), 1. / exp(ln_tau_O))(omega_s.col(0));

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

  // spatial effect
  array<Type> omega_s_A(n_i, 1);
  omega_s_A.setZero();
  omega_s_A.col(0) = A_st * vector<Type>(omega_s.col(0));

  vector<Type> pred(n_i);

  for (int i = 0; i < n_i; i++) {
    //eta_i(i,0) = omega_s_A(i,0) + B0;
    pred(i) = omega_s_A(i,0) + B0;
    jnll -= dnorm(y(i), pred(i), exp(ln_sigma), true);
  }
  // add prior penalty
  jnll -= penalty;

  // add jacobian adjustments
  jnll -= log(sqrt(8.0)) - log(pow(range,2.0)); // P(sigma)
  Type C = sqrt(exp(lgamma(1.0 + dhalf)) * pow(4*M_PI, dhalf));
  jnll -= log(C) + ln_kappa;

  // from calc_rf
  Type sigma_O = 1 / sqrt(Type(4.) * M_PI * exp(Type(2.) * ln_tau_O + Type(2.) * ln_kappa));
  REPORT(pred);
  ADREPORT(pred);
  REPORT(range);
  ADREPORT(range);
  REPORT(sigma_O);
  ADREPORT(sigma_O);
  return jnll;
}

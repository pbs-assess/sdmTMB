#define TMB_LIB_INIT R_init_sdmTMB
#include <TMB.hpp>

// ------------------ Main TMB template ----------------------------------------

template <class Type>
Type objective_function<Type>::operator()()
{
  using namespace R_inla;
  using namespace density;
  using namespace Eigen;

  // Vectors of real data
  DATA_VECTOR(y_i);      // response
  DATA_FACTOR(s_i);   // Random effect index for observation i
  DATA_INTEGER(n_s);  // number of sites (grids)
  DATA_STRUCT(spde, spde_t);

  // Parameters
  // Fixed effects

  PARAMETER(ln_tau_O);    // spatial process
  PARAMETER(ln_kappa);    // Matern parameter


  PARAMETER(ln_phi);           // sigma / dispersion / etc.

  PARAMETER_VECTOR(omega_s);    // spatial effects; n_s length

  // ------------------ End of parameters --------------------------------------

  int n_i = y_i.size();   // number of observations

  Type nll_data = 0;     // likelihood of data
  Type nll_omega = 0;    // spatial effects


  // ------------------ Geospatial ---------------------------------------------

  // Matern:
  Type range = sqrt(Type(8.0)) / exp(ln_kappa);


  Type sigma_O = 1 / sqrt(Type(4.0) * M_PI * exp(Type(2.0) * ln_tau_O) *
    exp(Type(2.0) * ln_kappa));
  REPORT(sigma_O);

  // Precision matrix
  Eigen::SparseMatrix<Type> Q;

  Q = R_inla::Q_spde(spde, exp(ln_kappa));

  // ------------------ Linear predictor----------------------------------------

  vector<Type> eta_i(n_i);
  for (int i = 0; i < n_i; i++) {
    eta_i(i) += omega_s(s_i(i));  // spatial
  }

  // ------------------ Probability of random effects --------------------------
  nll_omega += SCALE(GMRF(Q), 1.0 / exp(ln_tau_O))(omega_s);


  // ------------------ Probability of data given random effects ---------------

  // Type s1, s2;
  for (int i = 0; i < n_i; i++) {
    // if (!isNA(y_i(i))) {
    nll_data -= dnorm(y_i(i), eta_i(i), exp(ln_phi), true);
    // }
  }

  // ------------------ Projections --------------------------------------------

  // ---- Reporting ----------------------------------------------

  // REPORT(b_j)        // fixed effect parameters
  // REPORT(b_rw_t)     // fixed effect parameters
  REPORT(ln_tau_O);  // spatial process ln SD
  // REPORT(ln_tau_O_trend);  // spatial process ln SD
  // REPORT(ln_tau_E);  // spatio-temporal process ln SD
  // REPORT(ln_tau_V);  // spatio-temporal process ln SD
  // REPORT(sigma_E);
  REPORT(ln_phi);       // observation dispersion (depends on the distribution)
  // REPORT(thetaf);       // observation Tweedie mixing parameter
  // REPORT(epsilon_st);   // spatio-temporal effects; n_s by n_t matrix
  REPORT(omega_s);      // spatial effects; n_s length vector
  // REPORT(omega_s_trend);      // spatial effects; n_s length vector
  // REPORT(eta_fixed_i);  // fixed effect predictions in the link space
  // REPORT(eta_i);        // fixed and random effect predictions in link space
  REPORT(ln_kappa);     // Matern parameter
  REPORT(range);        // Matern approximate distance at 10% correlation

  // ------------------ Joint negative log likelihood --------------------------

  Type jnll = nll_data + nll_omega; // + nll_omega_trend + nll_varphi + nll_epsilon + nll_priors;
  return jnll;
}

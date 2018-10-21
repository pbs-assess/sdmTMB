#define TMB_LIB_INIT R_init_sdmTMB
#include <TMB.hpp>
#include "sdmTMB.h"

template <class Type>
bool isNA(Type x) {
  return R_IsNA(asDouble(x));
}

template <class Type>
vector<Type> Array2DToVector(array<Type> x) {
  x = x.transpose();
  int nr = x.rows();
  int nc = x.cols();
  vector<Type> res(nr * nc);
  for (int i = 0; i < nr; i++)
    for (int j = 0; j < nc; j++) {
      res[i * nc + j] = x(i, j);
    }
  return res;
}

template <class Type>
matrix<Type> MakeH(vector<Type> x) {
  matrix<Type> H(2, 2);
  H(0, 0) = exp(x(0));
  H(1, 0) = x(1);
  H(0, 1) = x(1);
  H(1, 1) = (1 + x(1) * x(1)) / exp(x(0));
  return H;
}

template <class Type>
vector<Type> Array1DToVector(array<Type> x) {
  int n = x.size();
  vector<Type> res(n);
  for (int i = 0; i < n; i++) res[i] = x(i);
  return res;
}

template <class Type>
vector<Type> RepeatVector(vector<Type> x, int times) {
  int n = x.size() * times;
  vector<Type> res(n);
  int k = 0;
  for (int i = 0; i < times; i++) {
    for (int j = 0; j < x.size(); j++) {
      res[k] = x(j);
      k++;
    }
  }
  return res;
}

enum valid_family {
  gaussian_family = 0,
  binomial_family = 1,
  tweedie_family  = 2,
  poisson_family  = 3,
  Gamma_family    = 4,
  nbinom2_family  = 5
};

enum valid_link {
  identity_link = 0,
  log_link      = 1,
  logit_link    = 2,
  inverse_link  = 3
};

template <class Type>
Type InverseLink(Type eta, int link) {
  Type out;
  switch (link) {
    case identity_link:
      out = eta;
      break;
    case log_link:
      out = exp(eta);
      break;
    case logit_link:
      out = eta;  // don't touch: we're using dbinom_robust() in logit space
      break; // FIXME make this more robust
    case inverse_link:
      out = Type(1.0) / eta;
      break;
    default:
      error("Link not implemented.");
  }
  return out;
}

// ------------------ Main TMB template ----------------------------------------

template <class Type>
Type objective_function<Type>::operator()() {
  using namespace R_inla;
  using namespace density;
  using namespace Eigen;

  // Vectors of real data
  DATA_VECTOR(y_i);   // response
  DATA_MATRIX(X_ij);  // model matrix

  DATA_FACTOR(s_i);   // Random effect index for observation i
  DATA_INTEGER(n_t);  // number of years
  DATA_INTEGER(n_s);  // number of sites (grids)

  // Indices for factors
  DATA_FACTOR(year_i);

  // Prediction?
  DATA_INTEGER(do_predict);

  // Distribution
  DATA_INTEGER(family);
  DATA_INTEGER(link);

  // SPDE objects from R-INLA
  DATA_STRUCT(spde_aniso, spde_aniso_t);
  DATA_STRUCT(spde, spde_t);
  PARAMETER_VECTOR(ln_H_input);
  DATA_INTEGER(anisotropy);

  // Projections
  DATA_SPARSE_MATRIX(proj_mesh);
  DATA_MATRIX(proj_X_ij);
  DATA_FACTOR(proj_year);

  // Spatial versus spatiotemporal
  DATA_INTEGER(spatial_only); //

  // ------------------ Parameters ---------------------------------------------

  // Parameters
  // Fixed effects
  PARAMETER_VECTOR(b_j);  // fixed effect parameters
  PARAMETER(ln_tau_O);    // spatial process
  PARAMETER(ln_tau_E);    // spatio-temporal process
  PARAMETER(ln_kappa);    // Matern parameter

  PARAMETER(thetaf);  // tweedie only
  PARAMETER(ln_phi);  // sigma / dispersion / etc.

  // Random effects
  // This is a matrix of spatial centers by years
  PARAMETER_VECTOR(omega_s);    // spatial effects; n_s length
  PARAMETER_ARRAY(epsilon_st);  // spatio-temporal effects; n_s by n_t matrix

  // ------------------ End of parameters --------------------------------------

  int n_i = y_i.size();  // number of observations
  int n_j = X_ij.cols();  // number of observations

  // Objective function is sum of negative log likelihood components
  Type nll_data = 0;  // likelihood of data
  Type nll_omega = 0;       // spatial effects
  Type nll_epsilon = 0;     // spatio-temporal effects
  Type nll_priors = 0;     // priors

  // ------------------ Priors -------------------------------------------------

  // nll_priors -= dnorm(ln_tau_O, Type(0.0), Type(10.0), true);
  // nll_priors -= dnorm(ln_tau_E, Type(0.0), Type(10.0), true);
  // nll_priors -= dnorm(ln_kappa, Type(0.0), Type(10.0), true);
  // nll_priors -= dnorm(ln_phi,   Type(0.0), Type(3.0), true);
  // for (int j = 0; j < n_j; j++)
  //   nll_priors -= dnorm(b_j(j), Type(0.0), Type(10.0), true);

  // ------------------ Geospatial ---------------------------------------------

  // Matern:
  Type range = sqrt(8.0) / exp(ln_kappa);
  Type sigma_O = 1 / sqrt(4 * M_PI * exp(2 * ln_tau_O) * exp(2 * ln_kappa));
  Type sigma_E = 1 / sqrt(4 * M_PI * exp(2 * ln_tau_E) * exp(2 * ln_kappa));

  // Precision matrix
  Eigen::SparseMatrix<Type> Q;
  if (anisotropy) {
    matrix<Type> H = MakeH(ln_H_input);
    Q = R_inla::Q_spde(spde_aniso, exp(ln_kappa), H);
    REPORT(H);
  }
  if (!anisotropy) {
    Q = R_inla::Q_spde(spde, exp(ln_kappa));
  }

  // ------------------ Linear predictor----------------------------------------

  vector<Type> eta_fixed_i = X_ij * b_j;
  vector<Type> mu_i(n_i), eta_i(n_i);
  for (int i = 0; i < n_i; i++) {
    eta_i(i) = eta_fixed_i(i) +                // fixed effects
               omega_s(s_i(i)) +               // spatial
               epsilon_st(s_i(i), year_i(i));  // spatio-temporal
    mu_i(i) = InverseLink(eta_i(i), link);
  }

  // ------------------ Probability of random effects --------------------------

  // Spatial effects:
  nll_omega += SCALE(GMRF(Q), 1.0 / exp(ln_tau_O))(omega_s);
  // Spatiotemporal effects:
  if (!spatial_only) {
    for (int t = 0; t < n_t; t++)
      nll_epsilon += SCALE(GMRF(Q), 1.0 / exp(ln_tau_E))(epsilon_st.col(t));
  }

  // ------------------ Probability of data given random effects ---------------

  Type s1, s2;
  for (int i = 0; i < n_i; i++) {
    if (!isNA(y_i(i))) {
      switch (family) {
        case gaussian_family:
          nll_data -= dnorm(y_i(i), mu_i(i), exp(ln_phi), true);
          break;
        case tweedie_family:
          s1 = invlogit(thetaf) + Type(1.0);
          nll_data -= dtweedie(y_i(i), mu_i(i), exp(ln_phi), s1, true);
          break;
        case binomial_family: // in logit space not inverse logit
          nll_data -= dbinom_robust(y_i(i), Type(1.0) /*size*/, mu_i(i), true);
          break;
        case poisson_family:
          nll_data -= dpois(y_i(i), mu_i(i), true);
          break;
        case Gamma_family:
          s1 = 1. / (pow(exp(ln_phi), 2.)); // s1=shape,ln_phi=CV,shape=1/CV^2
          nll_data -= dgamma(y_i(i), s1, mu_i(i) / s1, true);
          break;
        case nbinom2_family:
          error("Family not implemented.");
          // s1 = eta_i(i); // log(mu_i)
          // As in glmmTMB... FIXME honestly I'm not sure what's going on here:
          // s2 = 2. * s1 - log(pow(exp(ln_phi), 2.));     // log(var^2 - mu)
          // nll_data -= dnbinom_robust(y_i(i), s1, s2, true);
          break;
        default:
          error("Family not implemented.");
      }
    }
  }

  // ------------------ Projections --------------------------------------------

  if (do_predict) {
    vector<Type> proj_fe = proj_X_ij * b_j;
    vector<Type> proj_re_sp = proj_mesh * omega_s;
    vector<Type> proj_re_sp_st = RepeatVector(proj_re_sp, n_t);
    array<Type> proj_re_st(proj_mesh.rows(), n_t);
    for (int i = 0; i < n_t; i++)
      proj_re_st.col(i) = proj_mesh * Array1DToVector(epsilon_st.col(i));
    vector<Type> proj_re_st_vector = Array2DToVector(proj_re_st);
    vector<Type> proj_eta = proj_fe + proj_re_sp_st + proj_re_st_vector;
    REPORT(proj_fe);           // fixed effect projections
    REPORT(proj_re_sp);        // spatial random effect projections
    REPORT(proj_re_st_vector); // spatiotemporal random effect projections
    REPORT(proj_eta);          // combined projections (in link space)

  // ------------------ Derived quantities -------------------------------------
    vector<Type> total(n_t);
    for (int i = 0; i < proj_eta.size(); i++)  {
      total(proj_year(i)) += InverseLink(proj_eta(i), link);
    }
    vector<Type> log_total = log(total);
    REPORT(log_total);
    ADREPORT(log_total);
  }

  // ------------------ Reporting ----------------------------------------------

  REPORT(b_j)          // fixed effect parameters
  REPORT(ln_tau_O);    // spatial process ln SD
  REPORT(ln_tau_E);    // spatio-temporal process ln SD
  REPORT(sigma_O);
  REPORT(sigma_E);
  REPORT(ln_phi);      // observation dispersion (depends on the distribution)
  REPORT(thetaf);      // observation Tweedie mixing parameter
  REPORT(epsilon_st);  // spatio-temporal effects; n_s by n_t matrix
  REPORT(omega_s);     // spatial effects; n_s length vector
  REPORT(eta_fixed_i); // fixed effect predictions in the link space
  REPORT(eta_i);       // fixed and random effect predictions in link space
  REPORT(ln_kappa);    // Matern parameter
  REPORT(range);       // Matern approximate distance at 10% correlation

  // if (do_predict) ADREPORT(mu_i);

  // ------------------ Joint negative log likelihood --------------------------

  Type jnll = nll_data + nll_omega + nll_epsilon + nll_priors;
  return jnll;
}

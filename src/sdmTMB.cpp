#define TMB_LIB_INIT R_init_sdmTMB
#include <TMB.hpp>
#include "utils.h"
// #include <omp.h>

enum valid_family {
  gaussian_family = 0,
  binomial_family = 1,
  tweedie_family  = 2,
  poisson_family  = 3,
  Gamma_family    = 4,
  nbinom2_family  = 5,
  lognormal_family= 6,
  student_family  = 7,
  Beta_family     = 8,
  truncated_nbinom2_family  = 9,
  nbinom1_family  = 10,
  truncated_nbinom1_family  = 11
};

enum valid_link {
  identity_link = 0,
  log_link      = 1,
  logit_link    = 2,
  inverse_link  = 3
};

template <class Type>
Type InverseLink(Type eta, int link)
{
  Type out;
  switch (link) {
  case identity_link:
    out = eta;
    break;
  case log_link:
    out = exp(eta);
    break;
  case logit_link:
    out = invlogit(eta);
    break;
  case inverse_link:
    out = Type(1.0) / eta;
    break;
  default:
    error("Link not implemented.");
  }
  return out;
}

template <class Type>
Type Link(Type eta, int link)
{
  Type out;
  switch (link) {
  case identity_link:
    out = eta;
    break;
  case log_link:
    out = log(eta);
    break;
  case logit_link:
    out = logit(eta);
    break;
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
Type objective_function<Type>::operator()()
{
  using namespace R_inla;
  using namespace density;
  using namespace Eigen;

  // Set max number of OpenMP threads to help us optimize faster (as in glmmTMB)
  // max_parallel_regions = omp_get_max_threads();

  // Vectors of real data
  DATA_VECTOR(y_i);      // response
  DATA_MATRIX(X_ij);     // model matrix
  DATA_VECTOR(z_i);      // numeric vector for spatial covariate effect
  DATA_MATRIX(X_rw_ik);  // model matrix for random walk covariate(s)

  DATA_STRUCT(Zs, sdmTMB::LOM_t); // [L]ist [O]f (basis function matrices) [Matrices]
  DATA_STRUCT(proj_Zs, sdmTMB::LOM_t); // [L]ist [O]f (basis function matrices) [Matrices]
  DATA_MATRIX(Xs); // smoother linear effect matrix
  DATA_MATRIX(proj_Xs); // smoother linear effect matrix

  DATA_VECTOR_INDICATOR(keep, y_i); // https://rdrr.io/cran/TMB/man/oneStepPredict.html
  DATA_VECTOR(weights_i); // optional weights
  // DATA_VECTOR(offset_i); // optional offset

  DATA_INTEGER(n_t);  // number of years

  // Random intercepts:
  DATA_IMATRIX(RE_indexes);
  DATA_IMATRIX(proj_RE_indexes);
  DATA_IVECTOR(nobs_RE);
  DATA_IVECTOR(ln_tau_G_index);

  DATA_SPARSE_MATRIX(A); // INLA 'A' projection matrix for original data
  DATA_SPARSE_MATRIX(A_st); // INLA 'A' projection matrix for unique stations
  DATA_IVECTOR(A_spatial_index); // Vector of stations to match up A_st output

  // Indices for factors
  DATA_FACTOR(year_i);

  DATA_INTEGER(normalize_in_r);
  DATA_INTEGER(flag);
  DATA_INTEGER(share_range);

  // Prediction?
  DATA_INTEGER(do_predict);
  // With standard errors on the full projections?
  DATA_INTEGER(calc_se);
  // Should predictions be population (vs. individual-level) predictions?
  DATA_INTEGER(pop_pred);
  // Calculate total summed by year (e.g. biomass)?
  DATA_INTEGER(calc_time_totals);
  DATA_INTEGER(calc_quadratic_range);
  DATA_VECTOR(area_i); // area per prediction grid cell for index standardization

  DATA_VECTOR(priors_b_mean);
  DATA_MATRIX(priors_b_Sigma); // beta priors matrix
  DATA_INTEGER(priors_b_n);
  DATA_IVECTOR(priors_b_index);
  DATA_VECTOR(priors); // all other priors as a vector
  DATA_INTEGER(ar1_fields);
  DATA_INTEGER(rw_fields);
  DATA_INTEGER(include_spatial);
  DATA_INTEGER(random_walk);
  DATA_IVECTOR(exclude_RE);

  DATA_VECTOR(proj_lon);
  DATA_VECTOR(proj_lat);

  // Distribution
  DATA_INTEGER(family);
  DATA_INTEGER(link);
  DATA_SCALAR(df);  // Student-t DF
  DATA_VECTOR(size); // binomial, via glmmTMB

  // SPDE objects from R-INLA
  DATA_STRUCT(spde_aniso, spde_aniso_t);
  DATA_STRUCT(spde, spde_t);
  PARAMETER_VECTOR(ln_H_input);
  DATA_INTEGER(anisotropy);

  // Barrier
  DATA_INTEGER(barrier);
  DATA_STRUCT(spde_barrier, sdmTMB::spde_barrier_t);
  DATA_VECTOR(barrier_scaling); // scaling of range

  // Projections
  DATA_SPARSE_MATRIX(proj_mesh);
  DATA_MATRIX(proj_X_ij);
  DATA_MATRIX(proj_X_rw_ik);
  DATA_FACTOR(proj_year);
  DATA_VECTOR(proj_z_i);
  DATA_IVECTOR(proj_spatial_index);

  DATA_INTEGER(spatial_only);
  DATA_INTEGER(spatial_covariate);

  DATA_VECTOR(X_threshold);
  DATA_VECTOR(proj_X_threshold);
  DATA_INTEGER(threshold_func);
  // optional model for nonstationary st variance
  DATA_INTEGER(est_epsilon_model);
  DATA_VECTOR(epsilon_predictor);

  // optional stuff for penalized regression splines
  DATA_INTEGER(has_smooths);  // whether or not smooths are included
  DATA_IVECTOR(b_smooth_start);

  DATA_INTEGER(sim_re);
  // ------------------ Parameters ---------------------------------------------

  // Parameters
  // Fixed effects
  PARAMETER_VECTOR(b_j);  // fixed effect parameters
  PARAMETER_VECTOR(bs); // smoother linear effects
  PARAMETER(ln_tau_O);    // spatial process
  PARAMETER(ln_tau_Z);    // optional spatially varying covariate process
  PARAMETER(ln_tau_E);    // spatio-temporal process
  PARAMETER_VECTOR(ln_kappa);    // Matern parameter

  PARAMETER(thetaf);           // tweedie only
  PARAMETER(ln_phi);           // sigma / dispersion / etc.
  PARAMETER_VECTOR(ln_tau_V);  // random walk sigma
  PARAMETER(ar1_phi);          // AR1 fields correlation
  PARAMETER_VECTOR(ln_tau_G);  // random intercept sigmas
  PARAMETER_VECTOR(RE);        // random intercept deviations
  // Random effects
  PARAMETER_ARRAY(b_rw_t);  // random walk effects
  PARAMETER_VECTOR(omega_s);    // spatial effects; n_s length
  PARAMETER_VECTOR(zeta_s);    // spatial effects on covariate; n_s length
  PARAMETER_ARRAY(epsilon_st);  // spatio-temporal effects; n_s by n_t matrix
  PARAMETER_VECTOR(b_threshold);  // coefficients for threshold relationship (3)
  PARAMETER(b_epsilon); // slope coefficient for log-linear model on epsilon
  PARAMETER_VECTOR(b_smooth);  // P-spline smooth parameters
  PARAMETER_VECTOR(ln_smooth_sigma);  // variances of spline REs if included

  // Joint negative log-likelihood
  Type jnll = 0.0;

  // ------------------ End of parameters --------------------------------------

  int n_i = y_i.size();   // number of observations
  int n_RE = RE_indexes.cols();  // number of random effect intercepts

  // ------------------ Derived variables -------------------------------------------------
  Type s_slope, s_cut, s50, s95, s_max;
  // these are for linear model
  s_slope = b_threshold(0);
  s_cut = b_threshold(1);
  if (threshold_func == 2) {
    s50 = b_threshold(0); // threshold at which function is 50% of max
    s95 = b_threshold(0) + exp(b_threshold(1)); // threshold at which function is 95% of max
    s_max = b_threshold(2);
  }

  Type rho = sdmTMB::minus_one_to_one(ar1_phi);
  Type phi = exp(ln_phi);

  // ------------------ Geospatial ---------------------------------------------

  // Matern:
  vector<Type> range(2);
  range(0) = sqrt(Type(8.)) / exp(ln_kappa(0));
  range(1) = sqrt(Type(8.)) / exp(ln_kappa(1));

  if (include_spatial) {
    Type sigma_O = 1 / sqrt(Type(4.0) * M_PI * exp(Type(2.0) * ln_tau_O +
      Type(2.0) * ln_kappa(0)));
    Type sigma_Z = 1 / sqrt(Type(4.0) * M_PI * exp(Type(2.0) * ln_tau_Z +
      Type(2.0) * ln_kappa(0)));
    Type log_sigma_O = log(sigma_O);
    ADREPORT(log_sigma_O);
    REPORT(sigma_O);
    ADREPORT(sigma_O);
    Type log_sigma_Z = log(sigma_Z);
    ADREPORT(log_sigma_Z);
    REPORT(sigma_Z);
    ADREPORT(sigma_Z);
  }

  // optional non-stationary model on epsilon
  vector<Type> sigma_E(n_t);
  vector<Type> ln_tau_E_vec(n_t);
  //Type b_epsilon;
  if (!est_epsilon_model) { // constant model
    for (int i = 0; i < n_t; i++) {
      sigma_E(i) = 1 / sqrt(Type(4.0) * M_PI * exp(Type(2.0) * ln_tau_E + Type(2.0) * ln_kappa(1)));
      ln_tau_E_vec(i) = ln_tau_E;
    }
  }
  if (est_epsilon_model) { // loglinear model
    // epsilon_intcpt is the intercept parameter, derived from ln_tau_E.
    // For models with time as covariate, this is interpreted as sigma when covariate = 0.
    Type epsilon_intcpt = 1 / sqrt(Type(4.0) * M_PI * exp(Type(2.0) * ln_tau_E + Type(2.0) * ln_kappa(1)));
    Type log_epsilon_intcpt = log(epsilon_intcpt);
    Type log_epsilon_temp = 0.0;
    Type epsilon_cnst = - log(Type(4.0) * M_PI) / Type(2.0) - ln_kappa(1);
    for(int i = 0; i < n_t; i++) {
      log_epsilon_temp = log_epsilon_intcpt + b_epsilon * epsilon_predictor(i);
      sigma_E(i) = exp(log_epsilon_temp); // log-linear model
      ln_tau_E_vec(i) = -log_epsilon_temp + epsilon_cnst;
    }
  }

  Eigen::SparseMatrix<Type> Q_s; // Precision matrix
  Eigen::SparseMatrix<Type> Q_st; // Precision matrix

  if (barrier) {
    Q_s = Q_spde(spde_barrier, exp(ln_kappa(0)), barrier_scaling);
    if (!share_range) Q_st = Q_spde(spde_barrier, exp(ln_kappa(1)), barrier_scaling);
  } else {
    if (anisotropy) {
      matrix<Type> H = sdmTMB::MakeH(ln_H_input);
      Q_s = R_inla::Q_spde(spde_aniso, exp(ln_kappa(0)), H);
      if (!share_range) Q_st = R_inla::Q_spde(spde_aniso, exp(ln_kappa(1)), H);
      REPORT(H);
    }
    if (!anisotropy) {
      Q_s = R_inla::Q_spde(spde, exp(ln_kappa(0)));
      if (!share_range) Q_st = R_inla::Q_spde(spde, exp(ln_kappa(1)));
    }
  }
  if (share_range) Q_st = Q_s;

  // ------------------ Probability of random effects --------------------------

  // IID random intercepts:
  for (int g = 0; g < RE.size(); g++) {
    jnll -= dnorm(RE(g), Type(0.0), exp(ln_tau_G(ln_tau_G_index(g))), true);
    if (sim_re) SIMULATE{RE(g) = rnorm(Type(0), exp(ln_tau_G(ln_tau_G_index(g))));
    }
  }

  // Random walk effects (dynamic regression):
  if (random_walk) {
    for (int k = 0; k < X_rw_ik.cols(); k++) {
      // flat prior on the initial value... then:
      for (int t = 1; t < n_t; t++) {
        jnll += -dnorm(b_rw_t(t, k), b_rw_t(t - 1, k), exp(ln_tau_V(k)), true);
        if (sim_re) SIMULATE{b_rw_t(t, k) = rnorm(b_rw_t(t - 1, k), exp(ln_tau_V(k)));}
      }
    }
  }

  bool s = true;
  if (normalize_in_r) s = false;

  // Spatial (intercept) random effects:
  if (include_spatial) {
    // jnll += SCALE(GMRF(Q_s, s), 1. / exp(ln_tau_O))(omega_s);
    jnll += SCALE(GMRF(Q_s, s), 1. / exp(ln_tau_O))(omega_s);
    if (sim_re) {
      SIMULATE {
        GMRF(Q_s, s).simulate(omega_s);
        omega_s *= 1. / exp(ln_tau_O);
      }
    }
    if (spatial_covariate) {
      jnll += SCALE(GMRF(Q_s, s), 1. / exp(ln_tau_Z))(zeta_s);
      if (sim_re) {
        SIMULATE {
          GMRF(Q_s, s).simulate(zeta_s);
          zeta_s *= 1. / exp(ln_tau_Z);
        }
      }
    }
  }

  // Spatiotemporal random effects:
  if (!spatial_only) {
    if (!ar1_fields && !rw_fields) {
      for (int t = 0; t < n_t; t++)
        jnll += SCALE(GMRF(Q_st, s), 1. / exp(ln_tau_E_vec(t)))(epsilon_st.col(t));
      if (sim_re) {
        for (int t = 0; t < n_t; t++) { // untested!!
          vector<Type> epsilon_st_tmp(epsilon_st.rows());
          SIMULATE {GMRF(Q_st, s).simulate(epsilon_st_tmp);}
          epsilon_st.col(t) = epsilon_st_tmp / exp(ln_tau_E);
        }
      }
    } else {
      if (est_epsilon_model) { // time-varying epsilon sd
        if (sim_re) error("Simulation not implemented for time-varying epsilon SD yet.");
        if (ar1_fields) {
          jnll += SCALE(GMRF(Q_st, s), 1./exp(ln_tau_E_vec(0)))(epsilon_st.col(0));
          for (int t = 1; t < n_t; t++) {
            jnll += SCALE(GMRF(Q_st, s), 1./exp(ln_tau_E_vec(t)))((epsilon_st.col(t) -
              rho * epsilon_st.col(t - 1))/sqrt(1. - rho * rho));
          }
          int n_rows=epsilon_st.cols();
          int m_cols=epsilon_st.size()/n_rows;
          // This penalty added to match Kasper's AR1_t() implementation
          jnll += Type((n_rows-1)*m_cols) * log(sqrt(Type(1)-rho*rho));
        } else if (rw_fields) {
          jnll += SCALE(GMRF(Q_st, s), 1./exp(ln_tau_E_vec(0)))(epsilon_st.col(0));
          for (int t = 1; t < n_t; t++) {
            jnll += SCALE(GMRF(Q_st, s), 1./exp(ln_tau_E_vec(t)))(epsilon_st.col(t) - epsilon_st.col(t - 1));
          }
        } else {
          error("Field type not implemented.");
        }
      } else { // constant epsilon sd, keep calculations as is
        if (ar1_fields) {
          jnll += SCALE(SEPARABLE(AR1(rho), GMRF(Q_st, s)), 1./exp(ln_tau_E))(epsilon_st);
          if (sim_re) {
            SIMULATE {SEPARABLE(AR1(rho), GMRF(Q_st, s)).simulate(epsilon_st);}
            epsilon_st *= 1./exp(ln_tau_E);
          }
        } else if (rw_fields) {
          jnll += SCALE(GMRF(Q_st, s), 1./exp(ln_tau_E))(epsilon_st.col(0));
          for (int t = 1; t < n_t; t++) {
            jnll += SCALE(GMRF(Q_st, s), 1./exp(ln_tau_E))(epsilon_st.col(t) - epsilon_st.col(t - 1));
          }
          if (sim_re) {
            for (int t = 0; t < n_t; t++) { // untested!!
              vector<Type> epsilon_st_tmp(epsilon_st.rows());
              SIMULATE {GMRF(Q_st, s).simulate(epsilon_st_tmp);}
              epsilon_st_tmp *= 1./exp(ln_tau_E);
              if (t == 0) {
                epsilon_st.col(0) = epsilon_st_tmp;
              } else {
                epsilon_st.col(t) = epsilon_st.col(t-1) + epsilon_st_tmp;
              }
            }
          }
        } else {
          error("Field type not implemented.");
        }
      }
    }
  }
  if (flag == 0) return jnll;

  // ------------------ INLA projections ---------------------------------------

  // Here we are projecting the spatiotemporal and spatial random effects to the
  // locations of the data using the INLA 'A' matrices.
  array<Type> epsilon_st_A(A_st.rows(), n_t);
  for (int i = 0; i < n_t; i++)
    epsilon_st_A.col(i) = A_st * vector<Type>(epsilon_st.col(i));
  vector<Type> omega_s_A = A * omega_s;
  vector<Type> zeta_s_A = A * zeta_s;
  vector<Type> epsilon_st_A_vec(n_i);

  // ------------------ Linear predictor ---------------------------------------

  vector<Type> eta_fixed_i = X_ij * b_j;

  // p-splines/smoothers
  vector<Type> eta_smooth_i(X_ij.rows());
  eta_smooth_i.setZero();
  if (has_smooths) {
    for (int s = 0; s < b_smooth_start.size(); s++) { // iterate over # of smooth elements
      vector<Type> beta_s(Zs(s).cols());
      beta_s.setZero();
      for (int j = 0; j < beta_s.size(); j++) {
        beta_s(j) = b_smooth(b_smooth_start(s) + j);
        jnll -= dnorm(beta_s(j), Type(0), exp(ln_smooth_sigma(s)), true);
        if (sim_re) SIMULATE{beta_s(j) = rnorm(Type(0), exp(ln_smooth_sigma(s)));}
      }
      eta_smooth_i += Zs(s) * beta_s;
    }
    eta_smooth_i += Xs * bs;
  }

  // add threshold effect if specified
  if (threshold_func > 0) {
    if (threshold_func == 1) {
      // linear
      for (int i = 0; i < n_i; i++) {
        eta_fixed_i(i) += sdmTMB::linear_threshold(X_threshold(i), s_slope, s_cut);
      }
    } else {
      // logistic
      for (int i = 0; i < n_i; i++) {
        eta_fixed_i(i) += sdmTMB::logistic_threshold(X_threshold(i), s50, s95, s_max);
      }
    }
  }

  vector<Type> mu_i(n_i), eta_i(n_i), eta_rw_i(n_i), eta_iid_re_i(n_i);
  eta_rw_i.setZero();
  eta_iid_re_i.setZero();
  mu_i.setZero();
  eta_i.setZero();

  for (int i = 0; i < n_i; i++) {
    eta_i(i) = eta_fixed_i(i) + eta_smooth_i(i); // + offset_i(i);
    if (random_walk) {
      for (int k = 0; k < X_rw_ik.cols(); k++) {
        eta_rw_i(i) += X_rw_ik(i, k) * b_rw_t(year_i(i), k); // record it
        eta_i(i) += eta_rw_i(i);
      }
    }

    // Spatially varying effects:
    if (include_spatial) {
      eta_i(i) += omega_s_A(i);  // spatial
      if (spatial_covariate)
        eta_i(i) += zeta_s_A(i) * z_i(i); // spatially varying covariate
    }
    epsilon_st_A_vec(i) = epsilon_st_A(A_spatial_index(i), year_i(i)); // record it
    eta_i(i) += epsilon_st_A_vec(i); // spatiotemporal

    // IID random intercepts:
    int temp = 0;
    for (int k = 0; k < n_RE; k++) {
      if (k == 0) eta_iid_re_i(i) += RE(RE_indexes(i, k)); // record it
      if (k > 0) {
        temp += nobs_RE(k - 1);
        eta_iid_re_i(i) += RE(RE_indexes(i, k) + temp); // record it
      }
    }
    eta_i(i) += eta_iid_re_i(i);

    if (family == 1 && link == 2) {
      // binomial(link = "logit"); don't touch (using robust density function in logit space)
      mu_i(i) = eta_i(i);
    } else {
      mu_i(i) = InverseLink(eta_i(i), link);
    }
  }

  // ------------------ Probability of data given random effects ---------------

  // from glmmTMB:
  // close to zero: use for count data (cf binomial()$initialize)
#define zt_lik_nearzero(x,loglik_exp) ((x < Type(0.001)) ? -INFINITY : loglik_exp)

  Type s1, s2, s3, lognzprob, tmp_ll;
  REPORT(phi);
  ADREPORT(phi);
  for (int i = 0; i < n_i; i++) {
    if (!sdmTMB::isNA(y_i(i))) {
      switch (family) {
      case gaussian_family:
        tmp_ll = dnorm(y_i(i), mu_i(i), phi, true);
        SIMULATE{y_i(i) = rnorm(mu_i(i), phi);}
        break;
      case tweedie_family:
        s1 = invlogit(thetaf) + Type(1.0);
        if (!sdmTMB::isNA(priors(12))) jnll -= dnorm(s1, priors(12), priors(13), true);
        tmp_ll = dtweedie(y_i(i), mu_i(i), phi, s1, true);
        SIMULATE{y_i(i) = rtweedie(mu_i(i), phi, s1);}
        break;
      case binomial_family:  // in logit space not inverse logit
        tmp_ll = dbinom_robust(y_i(i), size(i), mu_i(i), true);
        SIMULATE{y_i(i) = rbinom(size(i), mu_i(i));}
        break;
      case poisson_family:
        tmp_ll = dpois(y_i(i), mu_i(i), true);
        SIMULATE{y_i(i) = rpois(mu_i(i));}
        break;
      case Gamma_family:
        s1 = exp(ln_phi);         // shape
        s2 = mu_i(i) / s1;        // scale
        tmp_ll = dgamma(y_i(i), s1, s2, true);
        SIMULATE{y_i(i) = rgamma(s1, s2);}
        // s1 = Type(1) / (pow(phi, Type(2)));  // s1=shape, ln_phi=CV,shape=1/CV^2
        // tmp_ll = dgamma(y_i(i), s1, mu_i(i) / s1, true);
        break;
      case nbinom2_family:
        s1 = log(mu_i(i)); // log(mu_i)
        s2 = 2. * s1 - ln_phi; // log(var - mu)
        tmp_ll = dnbinom_robust(y_i(i), s1, s2, true);
        SIMULATE { // from glmmTMB
          s1 = mu_i(i);
          s2 = mu_i(i) * (Type(1) + mu_i(i) / phi);
          y_i(i) = rnbinom2(s1, s2);
        }
        break;
      case truncated_nbinom2_family:
        s1 = log(mu_i(i)); // log(mu_i)
        s2 = 2. * s1 - ln_phi; // log(var - mu)
        tmp_ll = dnbinom_robust(y_i(i), s1, s2, true);
        s3 = logspace_add(Type(0), s1 - ln_phi);
        lognzprob = logspace_sub(Type(0), -phi * s3);
        tmp_ll -= lognzprob;
        tmp_ll = zt_lik_nearzero(y_i(i), tmp_ll); // from glmmTMB
        SIMULATE{y_i(i) = rtruncated_nbinom(asDouble(phi), 0, asDouble(mu_i(i)));}
        break;
      case nbinom1_family:
        s1 = log(mu_i(i));
        s2 = s1 + ln_phi;
        tmp_ll = dnbinom_robust(y_i(i), s1, s2, true);
        SIMULATE {y_i(i) = rnbinom2(mu_i(i), mu_i(i) * (Type(1) + phi));}
        break;
      case truncated_nbinom1_family:
        s1 = log(mu_i(i));
        s2 = s1 + ln_phi;
        tmp_ll = dnbinom_robust(y_i(i), s1, s2, true);
        s3 = logspace_add(Type(0), ln_phi);
        lognzprob = logspace_sub(Type(0), -mu_i(i) / phi * s3); // 1-prob(0)
        tmp_ll -= lognzprob;
        tmp_ll = zt_lik_nearzero(y_i(i), tmp_ll);
        SIMULATE{y_i(i) = rtruncated_nbinom(asDouble(mu_i(i)/phi), 0, asDouble(mu_i(i)));}
        break;
      case lognormal_family:
        tmp_ll = sdmTMB::dlnorm(y_i(i), log(mu_i(i)) - pow(phi, Type(2)) / Type(2), phi, true);
        SIMULATE{y_i(i) = exp(rnorm(log(mu_i(i)) - pow(phi, Type(2)) / Type(2), phi));}
        break;
      case student_family:
        tmp_ll = sdmTMB::dstudent(y_i(i), mu_i(i), exp(ln_phi), df, true);
        SIMULATE{y_i(i) = mu_i(i) + phi * rt(df);}
        break;
      case Beta_family: // Ferrari and Cribari-Neto 2004; betareg package
        s1 = mu_i(i) * phi;
        s2 = (Type(1) - mu_i(i)) * phi;
        tmp_ll = dbeta(y_i(i), s1, s2, true);
        SIMULATE{y_i(i) = rbeta(s1, s2);}
        break;
      default:
        error("Family not implemented.");
      }
      tmp_ll *= weights_i(i);
      jnll -= keep(i) * tmp_ll;
    }
  }

  // ------------------ Priors -------------------------------------------------

  // Construct special object for MVN distribution; always has mean 0.
  MVNORM_t<Type> neg_log_dmvnorm(priors_b_Sigma);
  // Apply nll on residual. Note that other univariate densities are positive
  // log-likelihoods but the dmvnorm is negative.
  // We're accumulating the neg LL, which is why this is a + sign.
  if(priors_b_n > 0) {
    vector<Type> b_j_subset(priors_b_n),b_mean_subset(priors_b_n);
    for(int j = 0; j < priors_b_n; j++) {
      b_j_subset(j) = b_j(priors_b_index(j));
      b_mean_subset(j) = priors_b_mean(j);
    }
    jnll += neg_log_dmvnorm(b_j_subset - b_mean_subset);
  }

  // start vector of priors:
  if (!sdmTMB::isNA(priors(0)) && !sdmTMB::isNA(priors(1)) && !sdmTMB::isNA(priors(2)) && !sdmTMB::isNA(priors(3))) {
    // std::cout << "Using spatial PC prior" << "\n";
    jnll -= sdmTMB::pc_prior_matern(ln_tau_O, ln_kappa(0), priors(0), priors(1), priors(2), priors(3), true);
  }
  if (!sdmTMB::isNA(priors(4)) && !sdmTMB::isNA(priors(5)) && !sdmTMB::isNA(priors(6)) && !sdmTMB::isNA(priors(7))) {
    // std::cout << "Using spatiotemporal PC prior" << "\n";
    jnll -= sdmTMB::pc_prior_matern(ln_tau_E, ln_kappa(1), priors(4), priors(5), priors(6), priors(7), true);
  }
  if (!sdmTMB::isNA(priors(8))) jnll -= dnorm(phi, priors(8), priors(9), true);
  if (!sdmTMB::isNA(priors(10))) jnll -= dnorm(rho, priors(10), priors(11), true);

  // Jacobians for Stan:
  // FIXME

  // ------------------ Predictions on new data --------------------------------

  if (do_predict) {
    vector<Type> proj_fe = proj_X_ij * b_j;
    // add threshold effect if specified
    if (threshold_func > 0) {
      if (threshold_func == 1) {
        // linear
        for (int i = 0; i < proj_X_ij.rows(); i++) {
          proj_fe(i) = proj_fe(i) + sdmTMB::linear_threshold(proj_X_threshold(i), s_slope, s_cut);
        }
      } else {
        // logistic
        for (int i = 0; i < proj_X_ij.rows(); i++) {
          proj_fe(i) = proj_fe(i) + sdmTMB::logistic_threshold(proj_X_threshold(i), s50, s95, s_max);
        }
      }
    }

    // Smoothers:
    vector<Type> proj_smooth_i(proj_X_ij.rows());
    proj_smooth_i.setZero();
    if (has_smooths) {
      for (int s = 0; s < b_smooth_start.size(); s++) { // iterate over # of smooth elements
        vector<Type> beta_s(proj_Zs(s).cols());
        beta_s.setZero();
        for (int j = 0; j < beta_s.size(); j++) {
          beta_s(j) = b_smooth(b_smooth_start(s) + j);
        }
        proj_smooth_i += proj_Zs(s) * beta_s;
      }
      proj_smooth_i += proj_Xs * bs;
    }
    for (int i = 0; i < proj_X_ij.rows(); i++) {
      proj_fe(i) += proj_smooth_i(i);
    }

    // IID random intercepts:
    vector<Type> proj_iid_re_i(proj_X_ij.rows());
    proj_iid_re_i.setZero();
    for (int i = 0; i < proj_X_ij.rows(); i++) {
      int temp = 0;
      for (int k = 0; k < n_RE; k++) {
        if (k == 0 && !exclude_RE(0)) proj_iid_re_i(i) += RE(proj_RE_indexes(i, k));
        if (k > 0) {
          temp += nobs_RE(k - 1);
          if (!exclude_RE(k)) proj_iid_re_i(i) += RE(proj_RE_indexes(i, k) + temp);
        }
      }
      proj_fe(i) += proj_iid_re_i(i);
    }

    // Random walk covariates:
    vector<Type> proj_rw_i(proj_X_ij.rows());
    proj_rw_i.setZero();
    if (random_walk) {
      for (int i = 0; i < proj_X_rw_ik.rows(); i++) {
        for (int k = 0; k < proj_X_rw_ik.cols(); k++) {
          proj_rw_i(i) += proj_X_rw_ik(i, k) * b_rw_t(proj_year(i), k);
          proj_fe(i) += proj_rw_i(i);
        }
      }
    }

    // Spatial and spatiotemporal random fields:
    vector<Type> proj_re_sp = proj_mesh * omega_s;
    vector<Type> proj_re_sp_st_all = sdmTMB::RepeatVector(proj_re_sp, n_t);
    array<Type> proj_re_st_temp(proj_mesh.rows(), n_t);
    array<Type> proj_re_st(proj_mesh.rows(), n_t);
    for (int i = 0; i < n_t; i++) {
      proj_re_st_temp.col(i) = proj_mesh * vector<Type>(epsilon_st.col(i));
      proj_re_st.col(i) = proj_re_st_temp.col(i);
    }

    // Spatially varying coefficients:
    vector<Type> proj_re_sp_cov(proj_X_ij.rows());
    vector<Type> proj_re_sp_slopes(proj_X_ij.rows());
    proj_re_sp_cov.setZero();
    proj_re_sp_slopes.setZero();
    if (spatial_covariate) {
      vector<Type> proj_re_sp_slopes_all = proj_mesh * zeta_s;
      for (int i = 0; i < proj_X_ij.rows(); i++) {
        proj_re_sp_cov(i) = proj_re_sp_slopes_all(proj_spatial_index(i)) * proj_z_i(i);
        proj_re_sp_slopes(i) = proj_re_sp_slopes_all(proj_spatial_index(i));
      }
    }

    // Pick out the appropriate spatial and/or or spatiotemporal values:
    vector<Type> proj_re_st_vector(proj_X_ij.rows());
    vector<Type> proj_re_sp_st(proj_X_ij.rows());
    proj_re_st_vector.setZero();
    proj_re_sp_st.setZero();
    for (int i = 0; i < proj_X_ij.rows(); i++) {
      proj_re_sp_st(i) = proj_re_sp_st_all(proj_spatial_index(i));
      proj_re_st_vector(i) = proj_re_st(proj_spatial_index(i), proj_year(i));
    }

    vector<Type> proj_rf = proj_re_sp_st + proj_re_st_vector + proj_re_sp_cov;
    vector<Type> proj_eta = proj_fe + proj_rf; // proj_fe includes RW and IID random effects

    REPORT(proj_fe);            // fixed effect projections
    REPORT(proj_re_sp_st);      // spatial random effect projections
    REPORT(proj_re_st_vector);  // spatiotemporal random effect projections
    REPORT(proj_re_sp_slopes);  // spatial slope projections
    REPORT(proj_re_sp_cov);     // spatial covariate projections
    REPORT(proj_eta);           // combined projections (in link space)
    REPORT(proj_rf);            // combined random field projections
    REPORT(proj_rw_i);          // random walk projections
    REPORT(proj_iid_re_i);      // IID random intercept projections

    if (calc_se) {
      if (pop_pred) {
        ADREPORT(proj_fe);
      } else {
        ADREPORT(proj_eta);
      }
    }

    if (calc_time_totals) {
      // ------------------ Derived quantities ---------------------------------

      // Total biomass etc.:
      vector<Type> total(n_t);
      for (int i = 0; i < proj_eta.size(); i++) {
        total(proj_year(i)) += InverseLink(proj_eta(i), link) * area_i(i);
      }
      vector<Type> link_total(n_t);
      for (int i = 0; i < n_t; i++) {
        link_total(i) = Link(total(i), link);
      }
      REPORT(link_total);
      ADREPORT(link_total);

      // Centre of gravity:
      vector<Type> cog_x(n_t);
      vector<Type> cog_y(n_t);
      for (int i = 0; i < proj_eta.size(); i++) {
        cog_x(proj_year(i)) += proj_lon(i) * InverseLink(proj_eta(i), link) * area_i(i);
        cog_y(proj_year(i)) += proj_lat(i) * InverseLink(proj_eta(i), link) * area_i(i);
      }
      for (int i = 0; i < n_t; i++) {
        cog_x(i) = cog_x(i) / total(i);
        cog_y(i) = cog_y(i) / total(i);
      }
      REPORT(cog_x);
      ADREPORT(cog_x);
      REPORT(cog_y);
      ADREPORT(cog_y);
    }
  }

  if (threshold_func == 1) { // linear breakpoint model
    REPORT(s_slope);
    ADREPORT(s_slope);
    REPORT(s_cut);
    ADREPORT(s_cut);
  }
  if (threshold_func == 2) { // logistic function model
    REPORT(s50);
    ADREPORT(s50);
    REPORT(s95);
    ADREPORT(s95);
    REPORT(s_max);
    ADREPORT(s_max);
  }
  if (calc_quadratic_range && b_j(1) < Type(0)) {
    vector<Type> quadratic_roots = sdmTMB::GetQuadraticRoots(b_j(1), b_j(0), Type(0.05));
    Type quadratic_low = quadratic_roots(0);
    Type quadratic_hi = quadratic_roots(1);
    Type quadratic_range = quadratic_roots(1) - quadratic_roots(0);
    if (quadratic_range < 0) quadratic_range = quadratic_range * -1.;
    Type quadratic_peak = quadratic_roots(2);
    Type quadratic_reduction = quadratic_roots(3);

    REPORT(quadratic_low);
    REPORT(quadratic_hi);
    REPORT(quadratic_range);
    REPORT(quadratic_peak);
    REPORT(quadratic_reduction);

    ADREPORT(quadratic_low);
    ADREPORT(quadratic_hi);
    ADREPORT(quadratic_range);
    ADREPORT(quadratic_peak);
    ADREPORT(quadratic_reduction);
  }
  if (est_epsilon_model) {
    REPORT(b_epsilon);
    ADREPORT(b_epsilon);
  }

  // ------------------ Reporting ----------------------------------------------

  vector<Type> log_sigma_E(n_t);
  for (int i = 0; i < n_t; i++) {
    log_sigma_E(i) = log(sigma_E(i));
  }
  ADREPORT(log_sigma_E);      // log spatio-temporal SD
  REPORT(sigma_E);      // spatio-temporal SD
  ADREPORT(sigma_E);      // spatio-temporal SD
  REPORT(epsilon_st_A_vec);   // spatio-temporal effects; vector
  REPORT(b_rw_t);   // time-varying effects
  REPORT(omega_s_A);      // spatial effects; n_s length vector
  REPORT(zeta_s_A);     // spatial covariate effects; n_s length vector
  REPORT(eta_fixed_i);  // fixed effect predictions in the link space
  REPORT(eta_smooth_i); // smooth effect predictions in the link space
  REPORT(eta_i);        // fixed and random effect predictions in link space
  REPORT(eta_rw_i);     // time-varying predictions in link space
  REPORT(eta_iid_re_i); // IID intercept random effect estimates
  REPORT(rho);          // AR1 correlation in -1 to 1 space
  REPORT(range);        // Matern approximate distance at 10% correlation
  ADREPORT(range);      // Matern approximate distance at 10% correlation
  vector<Type> log_range = log(range); // for SE
  ADREPORT(log_range);  // log Matern approximate distance at 10% correlation
  REPORT(b_smooth);     // smooth coefficients for penalized splines
  REPORT(ln_smooth_sigma); // standard deviations of smooth random effects, in log-space
  SIMULATE {
    REPORT(y_i);
    REPORT(omega_s)
    REPORT(omega_s_A)
    REPORT(epsilon_st)
    REPORT(epsilon_st_A_vec)
    REPORT(zeta_s)
    REPORT(zeta_s_A)
  }

  return jnll;
}

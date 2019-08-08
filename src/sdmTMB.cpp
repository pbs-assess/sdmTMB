#define TMB_LIB_INIT R_init_sdmTMB
#include <TMB.hpp>

template <class Type>
bool isNA(Type x)
{
  return R_IsNA(asDouble(x));
}

template <class Type>
Type dstudent(Type x, Type mean, Type sigma, Type df, int give_log = 0)
{
  // from metRology::dt.scaled()
  // dt((x - mean)/sd, df, ncp = ncp, log = TRUE) - log(sd)
  Type logres = dt((x - mean) / sigma, df, true) - log(sigma);
  if (give_log)
    return logres;
  else
    return exp(logres);
}

template <class Type>
Type dlnorm(Type x, Type meanlog, Type sdlog, int give_log = 0)
{
  Type logres = dnorm(log(x), meanlog, sdlog, true) - log(x);
  if (give_log)
    return logres;
  else
    return exp(logres);
}

template <class Type>
Type minus_one_to_one(Type x)
{
  return Type(2) * invlogit(x) - Type(1);
}

template <class Type>
matrix<Type> MakeH(vector<Type> x)
{
  matrix<Type> H(2, 2);
  H(0, 0) = exp(x(0));
  H(1, 0) = x(1);
  H(0, 1) = x(1);
  H(1, 1) = (1 + x(1) * x(1)) / exp(x(0));
  return H;
}

template <class Type>
vector<Type> RepeatVector(vector<Type> x, int times)
{
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

template <class Type>
vector<Type> GetQuadraticRoots(Type a, Type b, Type threshold)
{
  vector<Type> res(2);
  Type c = 1.; // doesn't matter; setting to an arbitrary value
  Type crit_y = (a * pow(-b / (2. * a), 2.) + b * (-b / (2. * a)) + c) + log(threshold);
  // solve for 0 = ax2 + bx + (c - crit_y)
  c = c - crit_y;
  res(0) = (b - sqrt(pow(b, 2.) - 4. * c * a))/(2.*a);
  res(1) = (b + sqrt(pow(b, 2.) - 4. * c * a))/(2.*a);
  return res;
}

enum valid_family {
  gaussian_family = 0,
  binomial_family = 1,
  tweedie_family  = 2,
  poisson_family  = 3,
  Gamma_family    = 4,
  nbinom2_family  = 5,
  lognormal_family= 6,
  student_family  = 7,
  Beta_family     = 8
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

// ------------------ Main TMB template ----------------------------------------

template <class Type>
Type objective_function<Type>::operator()()
{
  using namespace R_inla;
  using namespace density;
  using namespace Eigen;

  // Vectors of real data
  DATA_VECTOR(y_i);      // response
  DATA_MATRIX(X_ij);     // model matrix
  DATA_VECTOR(t_i);      // numeric year vector -- only for spatial_trend==1
  DATA_MATRIX(X_rw_ik);  // model matrix for random walk covariate(s)

  DATA_INTEGER(n_t);  // number of years

  DATA_SPARSE_MATRIX(A); // INLA 'A' projection matrix for original data
  DATA_SPARSE_MATRIX(A_st); // INLA 'A' projection matrix for unique stations
  DATA_IVECTOR(A_spatial_index); // Vector of stations to match up A_st output

  // Indices for factors
  DATA_FACTOR(year_i);

  DATA_INTEGER(flag); // flag=0 => only prior returned; used when normalizing in R
  DATA_INTEGER(normalize_in_r);

  // Prediction?
  DATA_INTEGER(do_predict);
  // With standard errors on the full projections?
  DATA_INTEGER(calc_se);
  // Calculate total summed by year (e.g. biomass)?
  DATA_INTEGER(calc_time_totals);
  DATA_INTEGER(calc_quadratic_range);

  DATA_INTEGER(enable_priors);
  DATA_INTEGER(ar1_fields);
  // DATA_INTEGER(separable_ar1);
  DATA_INTEGER(include_spatial);
  DATA_INTEGER(random_walk);

  DATA_VECTOR(proj_lon);
  DATA_VECTOR(proj_lat);

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
  DATA_MATRIX(proj_X_rw_ik);
  DATA_FACTOR(proj_year);
  DATA_VECTOR(proj_t_i);
  DATA_IVECTOR(proj_spatial_index);

  DATA_INTEGER(spatial_only);
  DATA_INTEGER(spatial_trend);

  // ------------------ Parameters ---------------------------------------------

  // Parameters
  // Fixed effects
  PARAMETER_VECTOR(b_j);  // fixed effect parameters
  PARAMETER(ln_tau_O);    // spatial process
  PARAMETER(ln_tau_O_trend);    // optional spatial process on the trend
  PARAMETER(ln_tau_E);    // spatio-temporal process
  PARAMETER(ln_kappa);    // Matern parameter

  PARAMETER(thetaf);           // tweedie only
  PARAMETER(ln_phi);           // sigma / dispersion / etc.
  PARAMETER_VECTOR(ln_tau_V);  // random walk sigma
  PARAMETER(ar1_phi);          // AR1 fields correlation

  // Random effects
  PARAMETER_ARRAY(b_rw_t);  // random walk effects
  PARAMETER_VECTOR(omega_s);    // spatial effects; n_s length
  PARAMETER_VECTOR(omega_s_trend);    // spatial effects on trend; n_s length
  PARAMETER_ARRAY(epsilon_st);  // spatio-temporal effects; n_s by n_t matrix

  // ------------------ End of parameters --------------------------------------

  int n_i = y_i.size();   // number of observations
  int n_j = X_ij.cols();  // number of observations

  Type nll_data = 0;     // likelihood of data
  Type nll_varphi = 0;   // random walk effects
  Type nll_omega = 0;    // spatial effects
  Type nll_omega_trend = 0;    // spatial trend effects
  Type nll_epsilon = 0;  // spatio-temporal effects
  Type nll_priors = 0;   // priors

  // ------------------ Priors -------------------------------------------------

  if (enable_priors) {
    nll_priors -= dnorm(ln_tau_O, Type(0.0), Type(1.0), true);
    nll_priors -= dnorm(ln_tau_E, Type(0.0), Type(1.0), true);
    nll_priors -= dnorm(ln_kappa, Type(0.0), Type(2.0), true);
    nll_priors -= dnorm(ln_phi, Type(0.0), Type(1.0), true);
    for (int j = 0; j < n_j; j++)
      nll_priors -= dnorm(b_j(j), Type(0.0), Type(5.0), true);
    if (spatial_trend) {
      nll_priors -= dnorm(ln_tau_O_trend, Type(0.0), Type(1.0), true);
    }
  }

  // ------------------ Geospatial ---------------------------------------------

  // Matern:
  Type range = sqrt(Type(8.0)) / exp(ln_kappa);

  if (include_spatial) {
    Type sigma_O = 1 / sqrt(Type(4.0) * M_PI * exp(Type(2.0) * ln_tau_O) *
                            exp(Type(2.0) * ln_kappa));
    Type sigma_O_trend = 1 / sqrt(Type(4.0) * M_PI * exp(Type(2.0) * ln_tau_O_trend) *
      exp(Type(2.0) * ln_kappa));
    REPORT(sigma_O);
    ADREPORT(sigma_O);
    REPORT(sigma_O_trend);
    ADREPORT(sigma_O_trend);
  }
  Type sigma_E = 1 / sqrt(Type(4.0) * M_PI * exp(Type(2.0) * ln_tau_E) *
                          exp(Type(2.0) * ln_kappa));

  Eigen::SparseMatrix<Type> Q; // Precision matrix
  if (anisotropy) {
    matrix<Type> H = MakeH(ln_H_input);
    Q = R_inla::Q_spde(spde_aniso, exp(ln_kappa), H);
    REPORT(H);
  }
  if (!anisotropy) {
    Q = R_inla::Q_spde(spde, exp(ln_kappa));
  }

  // ------------------ INLA projections ---------------------------------------

  // Here we are projecting the spatiotemporal and spatial random effects to the
  // locations of the data using the INLA 'A' matrices.
  array<Type> epsilon_st_A(A_st.rows(), n_t);
  for (int i = 0; i < n_t; i++)
    epsilon_st_A.col(i) = A_st * vector<Type>(epsilon_st.col(i));
  vector<Type> omega_s_A = A * omega_s;
  vector<Type> omega_s_trend_A = A * omega_s_trend;
  vector<Type> epsilon_st_A_vec(n_i);

  // ------------------ Linear predictor ---------------------------------------

  vector<Type> eta_fixed_i = X_ij * b_j;
  vector<Type> mu_i(n_i), eta_i(n_i);
  vector<Type> eta_rw_i(n_i);
  for (int i = 0; i < n_i; i++) {
    eta_i(i) = Type(0);
    eta_rw_i(i) = Type(0);
  }
  for (int i = 0; i < n_i; i++) {
    eta_i(i) = eta_fixed_i(i);
    if (random_walk)
      for (int k = 0; k < X_rw_ik.cols(); k++) {
        eta_rw_i(i) += X_rw_ik(i, k) * b_rw_t(year_i(i), k); // record it
        eta_i(i) += eta_rw_i(i);
      }
      if (include_spatial) {
        eta_i(i) += omega_s_A(i);  // spatial
        if (spatial_trend)
          eta_i(i) += omega_s_trend_A(i) * t_i(i); // spatial trend
      }
      epsilon_st_A_vec(i) = epsilon_st_A(A_spatial_index(i), year_i(i)); // record it
      eta_i(i) += epsilon_st_A_vec(i); // spatiotemporal

      if (family == 1 && link == 2) {
        // binomial(link = "logit"); don't touch (using robust density function in logit space)
        mu_i(i) = eta_i(i);
      } else {
        mu_i(i) = InverseLink(eta_i(i), link);
      }
  }

  // ------------------ Probability of random effects --------------------------

  // Random walk effects (dynamic regression):
  if (random_walk) {
    for (int k = 0; k < X_rw_ik.cols(); k++) {
      // flat prior on the initial value... then:
      for (int t = 1; t < n_t; t++) {
        nll_varphi += -dnorm(b_rw_t(t, k), b_rw_t(t - 1, k), exp(ln_tau_V(k)), true);
      }
    }
  }

  Type rho = minus_one_to_one(ar1_phi);
  bool s = true;
  if (normalize_in_r) s = false;

  // Spatial (intercept) random effects:
  if (include_spatial)
    nll_omega += SCALE(GMRF(Q, s), 1.0 / exp(ln_tau_O))(omega_s);
  // Spatial trend random effects:
  if (spatial_trend)
    nll_omega_trend += SCALE(GMRF(Q, s), 1.0 / exp(ln_tau_O_trend))(omega_s_trend);
  // Spatiotemporal random effects:
  if (!spatial_only) {
    if (!ar1_fields) {
      for (int t = 0; t < n_t; t++)
        nll_epsilon += SCALE(GMRF(Q, s), 1. / exp(ln_tau_E))(epsilon_st.col(t));
    } else {
      // if (!separable_ar1) {
      //   nll_epsilon += SCALE(GMRF(Q, s), 1./exp(ln_tau_E))(epsilon_st.col(0));
      //   for (int t = 1; t < n_t; t++) {
      //     nll_epsilon += SCALE(GMRF(Q, s), 1./exp(ln_tau_E))(epsilon_st.col(t) -
      //       rho * epsilon_st.col(t - 1));
      //   }
      // } else {
        nll_epsilon += SCALE(SEPARABLE(AR1(rho), GMRF(Q, s)), 1./exp(ln_tau_E))(epsilon_st);
      // }
    }
  }

  // Normalization of GMRFs during outer-optimization step in R:
  if (normalize_in_r) {
    Type nll_gmrf = nll_epsilon + nll_omega + nll_omega_trend;
    if (flag == 0) return(nll_gmrf);
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
        case binomial_family:  // in logit space not inverse logit
          nll_data -= dbinom_robust(y_i(i), Type(1.0) /*size*/, mu_i(i), true);
          break;
        case poisson_family:
          nll_data -= dpois(y_i(i), mu_i(i), true);
          break;
        case Gamma_family:
          s1 = Type(1) / (pow(exp(ln_phi), Type(2)));  // s1=shape,ln_phi=CV,shape=1/CV^2
          nll_data -= dgamma(y_i(i), s1, mu_i(i) / s1, true);
          break;
        case nbinom2_family:
          s1 = log(mu_i(i)); // log(mu_i)
          s2 = 2. * s1 - ln_phi; // log(var - mu)
          nll_data -= dnbinom_robust(y_i(i), s1, s2, true);
          break;
        case lognormal_family:
          nll_data -= dlnorm(y_i(i), mu_i(i) - pow(exp(ln_phi), Type(2)) / Type(2), exp(ln_phi), true);
          break;
        case student_family:
          nll_data -= dstudent(y_i(i), mu_i(i), exp(ln_phi), Type(3) /*df*/, true);
          break;
        case Beta_family: // Ferrari and Cribari-Neto 2004; betareg package
          s1 = mu_i(i) * exp(ln_phi);
          s2 = (Type(1) - mu_i(i)) * exp(ln_phi);
          nll_data -= dbeta(y_i(i), s1, s2, true);
          break;
        default:
          error("Family not implemented.");
      }
    }
  }

  // ------------------ Predictions on new data --------------------------------

  if (do_predict) {
    vector<Type> proj_fe = proj_X_ij * b_j;
    if (random_walk) {
      for (int i = 0; i < proj_X_rw_ik.rows(); i++) {
        for (int k = 0; k < proj_X_rw_ik.cols(); k++) {
          proj_fe(i) += proj_X_rw_ik(i, k) * b_rw_t(proj_year(i), k);
        }
      }
    }
    vector<Type> proj_re_sp = proj_mesh * omega_s;
    vector<Type> proj_re_sp_st_all = RepeatVector(proj_re_sp, n_t);
    array<Type> proj_re_st_temp(proj_mesh.rows(), n_t);
    array<Type> proj_re_st(proj_mesh.rows(), n_t);
    for (int i = 0; i < n_t; i++) {
      proj_re_st_temp.col(i) = proj_mesh * vector<Type>(epsilon_st.col(i));
      proj_re_st.col(i) = proj_re_st_temp.col(i);
    }

    vector<Type> proj_re_sp_trend(proj_X_ij.rows());
    vector<Type> proj_re_sp_slopes(proj_X_ij.rows());
    for (int i = 0; i < proj_X_ij.rows(); i++) {
      proj_re_sp_trend(i) = Type(0);
      proj_re_sp_slopes(i) = Type(0);
    }

    if (spatial_trend) {
      vector<Type> proj_re_sp_slopes_all = proj_mesh * omega_s_trend;
      for (int i = 0; i < proj_X_ij.rows(); i++) {
        proj_re_sp_trend(i) = proj_re_sp_slopes_all(proj_spatial_index(i)) * proj_t_i(i);
        proj_re_sp_slopes(i) = proj_re_sp_slopes_all(proj_spatial_index(i));
      }
    }

    // Pick out the appropriate spatial and/or or spatiotemporal values:
    vector<Type> proj_re_st_vector(proj_X_ij.rows());
    vector<Type> proj_re_sp_st(proj_X_ij.rows());
    for (int i = 0; i < proj_X_ij.rows(); i++) {
      proj_re_st_vector(i) = Type(0);
      proj_re_sp_st(i) = Type(0);
    }
    for (int i = 0; i < proj_X_ij.rows(); i++) {
      proj_re_sp_st(i) = proj_re_sp_st_all(proj_spatial_index(i));
      proj_re_st_vector(i) = proj_re_st(proj_spatial_index(i), proj_year(i));
    }

    vector<Type> proj_eta = proj_fe + proj_re_sp_st +
      proj_re_st_vector + proj_re_sp_trend;
    vector<Type> proj_rf = proj_re_sp_st + proj_re_st_vector + proj_re_sp_trend;
    REPORT(proj_fe);            // fixed effect projections
    REPORT(proj_re_sp_st);      // spatial random effect projections
    REPORT(proj_re_st_vector);  // spatiotemporal random effect projections
    REPORT(proj_re_sp_slopes);  // spatial slope projections
    REPORT(proj_re_sp_trend);   // spatial trend projections (slope * time)
    REPORT(proj_eta);           // combined projections (in link space)
    REPORT(proj_rf);            // combined random field projections

    if (calc_se) ADREPORT(proj_eta);

    if (calc_time_totals) {
      // ------------------ Derived quantities ---------------------------------

      // Total biomass:
      vector<Type> total(n_t);
      for (int i = 0; i < proj_eta.size(); i++) {
        total(proj_year(i)) += InverseLink(proj_eta(i), link);
      }
      vector<Type> log_total = log(total);
      REPORT(log_total);
      ADREPORT(log_total);

      // Centre of gravity:
      vector<Type> cog_x(n_t);
      vector<Type> cog_y(n_t);
      for (int i = 0; i < proj_eta.size(); i++) {
        cog_x(proj_year(i)) += proj_lon(i) * InverseLink(proj_eta(i), link);
        cog_y(proj_year(i)) += proj_lat(i) * InverseLink(proj_eta(i), link);
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

  if (calc_quadratic_range && b_j(1) < Type(0)) {
    vector<Type> quadratic_roots = GetQuadraticRoots(b_j(0), b_j(1), Type(0.05));
    Type quadratic_range = quadratic_roots(1) - quadratic_roots(0);
    if (quadratic_range < 0) quadratic_range = quadratic_range * -1.;
    REPORT(quadratic_roots);
    REPORT(quadratic_range);
    ADREPORT(quadratic_roots);
    ADREPORT(quadratic_range);
  }

  // ------------------ Reporting ----------------------------------------------

  REPORT(sigma_E);      // spatio-temporal process parameter
  ADREPORT(sigma_E);      // spatio-temporal process parameter
  REPORT(epsilon_st_A_vec);   // spatio-temporal effects; vector
  REPORT(b_rw_t);   // time-varying effects
  REPORT(omega_s_A);      // spatial effects; n_s length vector
  REPORT(omega_s_trend_A); // spatial trend effects; n_s length vector
  REPORT(eta_fixed_i);  // fixed effect predictions in the link space
  REPORT(eta_i);        // fixed and random effect predictions in link space
  REPORT(eta_rw_i);     // time-varying predictions in link space
  REPORT(rho);          // AR1 correlation in -1 to 1 space
  REPORT(range);        // Matern approximate distance at 10% correlation
  ADREPORT(range);      // Matern approximate distance at 10% correlation

  // ------------------ Joint negative log likelihood --------------------------

  Type jnll = nll_data + nll_omega + nll_omega_trend + nll_varphi + nll_epsilon + nll_priors;
  return jnll;
}

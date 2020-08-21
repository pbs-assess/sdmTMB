#include <TMB.hpp>
// #include <omp.h>

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
  vector<Type> res(4);
  Type c = 1.; // doesn't matter; setting to an arbitrary value
  Type crit_y = (a * pow(-b / (2. * a), 2.) + b * (-b / (2. * a)) + c) + log(threshold);
  // solve for 0 = ax2 + bx + (c - crit_y)
  c = c - crit_y;
  res(0) = -1. * (b - sqrt(pow(b, 2.) - 4. * c * a))/(2.*a);
  res(1) = -1. * (b + sqrt(pow(b, 2.) - 4. * c * a))/(2.*a);

  // calculate vertex of parabola
  Type xpeak = -b / (2.*a);
  // res(2) is the hi/lowpoint of parabola evaluated at xpeak
  res(2) = (a * (pow(xpeak, 2.)) + b * (xpeak) + c);

  // calculate reduction of changing from mean to +/- 1 SD
  res(3) = (a * (pow(xpeak+1, 2.)) + b * (xpeak+1) + c) / res(2);
  return res;
}

template <class Type>
Type linear_threshold(Type x, Type slope, Type cutpoint)
{
  // linear threshold model. relationship linear up to a point then constant
  // keep all parameters unconstrained - slope and scale can be neg/pos,
  // as can cutpoint if covariate is scaled ~ N(0,1).
  Type pred;
  if(x < cutpoint) {
    pred = x * slope;
  } else {
    pred = x * cutpoint;
  }
  return pred;
}

template <class Type>
Type logistic_threshold(Type x, Type s50, Type s95, Type scale)
{
  // logistic threshold model. similar to length or size based selectvitiy
  // in fisheries, parameterized by the points at which f(x) = 0.5, or 0.95
  // s50 and scale are unconstrained. s95 has to be > s50 though, so modeled as
  // s95 = s50 + exp(b(1))
  //Type s95 = s50 + exp(soffset); // this done outside function
  Type pred = (scale) * Type(1.0)/(Type(1.0) + exp(-log(Type(19.0)) * (x - s50) / (s95 - s50)));
  return pred;
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

  // Set max number of OpenMP threads to help us optimize faster (as in glmmTMB)
  // max_parallel_regions = omp_get_max_threads();

  // Vectors of real data
  DATA_VECTOR(y_i);      // response
  DATA_MATRIX(X_ij);     // model matrix
  DATA_VECTOR(t_i);      // numeric year vector -- only for spatial_trend==1
  DATA_MATRIX(X_rw_ik);  // model matrix for random walk covariate(s)

  DATA_VECTOR_INDICATOR(keep, y_i); // https://rdrr.io/cran/TMB/man/oneStepPredict.html
  DATA_VECTOR(weights_i); // optional weights
  // DATA_VECTOR(offset_i); // optional offset

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
  DATA_VECTOR(area_i); // area per prediction grid cell for index standardization

  DATA_INTEGER(enable_priors);
  DATA_INTEGER(ar1_fields);
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

  DATA_VECTOR(X_threshold);
  DATA_VECTOR(proj_X_threshold);
  DATA_INTEGER(threshold_func);
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

  PARAMETER_VECTOR(b_threshold);  // coefficients for threshold relationship (3)

  // Joint negative log-likelihood
  Type jnll = 0;

  // ------------------ End of parameters --------------------------------------

  int n_i = y_i.size();   // number of observations
  int n_j = X_ij.cols();  // number of observations

  // Type nll_data = 0;     // likelihood of data
  // Type nll_varphi = 0;   // random walk effects
  // Type nll_omega = 0;    // spatial effects
  // Type nll_omega_trend = 0;    // spatial trend effects
  // Type nll_epsilon = 0;  // spatio-temporal effects
  // Type nll_priors = 0;   // priors

  // ------------------ Derived variables -------------------------------------------------
  Type s_slope, s_cut, s50, s95, s_max;
  // these are for linear model
  s_slope = b_threshold(0);
  s_cut = b_threshold(1);
  if(threshold_func == 2) {
    s50 = b_threshold(0); // threshold at which function is 50% of max
    s95 = b_threshold(0) + exp(b_threshold(1)); // threshold at which function is 95% of max
    s_max = b_threshold(2);
  }
  // ------------------ Priors -------------------------------------------------

  if (enable_priors) {
    jnll -= dnorm(ln_tau_O, Type(0.0), Type(1.0), true);
    jnll -= dnorm(ln_tau_E, Type(0.0), Type(1.0), true);
    jnll -= dnorm(ln_kappa, Type(0.0), Type(2.0), true);
    jnll -= dnorm(ln_phi, Type(0.0), Type(1.0), true);
    for (int j = 0; j < n_j; j++)
      jnll -= dnorm(b_j(j), Type(0.0), Type(5.0), true);
    if (spatial_trend) {
      jnll -= dnorm(ln_tau_O_trend, Type(0.0), Type(1.0), true);
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
  // add threshold effect if specified
  if (threshold_func > 0) {
    if (threshold_func == 1) {
      // linear
      for (int i = 0; i < n_i; i++) {
        eta_fixed_i(i) += linear_threshold(X_threshold(i), s_slope, s_cut);
      }
    } else {
      // logistic
      for (int i = 0; i < n_i; i++) {
        eta_fixed_i(i) += logistic_threshold(X_threshold(i), s50, s95, s_max);
      }
    }
  }

  vector<Type> mu_i(n_i), eta_i(n_i);
  vector<Type> eta_rw_i(n_i);
  for (int i = 0; i < n_i; i++) {
    eta_i(i) = Type(0);
    eta_rw_i(i) = Type(0);
  }
  for (int i = 0; i < n_i; i++) {
    eta_i(i) = eta_fixed_i(i); // + offset_i(i);
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
        jnll += -dnorm(b_rw_t(t, k), b_rw_t(t - 1, k), exp(ln_tau_V(k)), true);
      }
    }
  }

  Type rho = minus_one_to_one(ar1_phi);
  bool s = true;
  if (normalize_in_r) s = false;

  // Spatial (intercept) random effects:
  if (include_spatial)
    jnll += SCALE(GMRF(Q, s), 1.0 / exp(ln_tau_O))(omega_s);
  // Spatial trend random effects:
  if (spatial_trend)
    jnll += SCALE(GMRF(Q, s), 1.0 / exp(ln_tau_O_trend))(omega_s_trend);
  // Spatiotemporal random effects:
  if (!spatial_only) {
    if (!ar1_fields) {
      for (int t = 0; t < n_t; t++)
        jnll += SCALE(GMRF(Q, s), 1. / exp(ln_tau_E))(epsilon_st.col(t));
    } else {
      // if (!separable_ar1) {
      //   nll_epsilon += SCALE(GMRF(Q, s), 1./exp(ln_tau_E))(epsilon_st.col(0));
      //   for (int t = 1; t < n_t; t++) {
      //     nll_epsilon += SCALE(GMRF(Q, s), 1./exp(ln_tau_E))(epsilon_st.col(t) -
      //       rho * epsilon_st.col(t - 1));
      //   }
      // } else {
      jnll += SCALE(SEPARABLE(AR1(rho), GMRF(Q, s)), 1./exp(ln_tau_E))(epsilon_st);
      // }
    }
  }

  // Normalization of GMRFs during outer-optimization step in R:
  if (normalize_in_r) {
    // Type nll_gmrf = nll_epsilon + nll_omega + nll_omega_trend;
    // if (flag == 0) return(jnll);
  }

  // ------------------ Probability of data given random effects ---------------

  Type s1, s2;
  for (int i = 0; i < n_i; i++) {
    if (!isNA(y_i(i))) {
      switch (family) {
        case gaussian_family:
          jnll -= keep(i) * dnorm(y_i(i), mu_i(i), exp(ln_phi), true) * weights_i(i);
          break;
        case tweedie_family:
          s1 = invlogit(thetaf) + Type(1.0);
          jnll -= keep(i) * dtweedie(y_i(i), mu_i(i), exp(ln_phi), s1, true) * weights_i(i);
          break;
        case binomial_family:  // in logit space not inverse logit
          jnll -= keep(i) * dbinom_robust(y_i(i), Type(1.0) /*size*/, mu_i(i), true) * weights_i(i);
          break;
        case poisson_family:
          jnll -= keep(i) * dpois(y_i(i), mu_i(i), true) * weights_i(i);
          break;
        case Gamma_family:
          s1 = Type(1) / (pow(exp(ln_phi), Type(2)));  // s1=shape,ln_phi=CV,shape=1/CV^2
          jnll -= keep(i) * dgamma(y_i(i), s1, mu_i(i) / s1, true) * weights_i(i);
          break;
        case nbinom2_family:
          s1 = log(mu_i(i)); // log(mu_i)
          s2 = 2. * s1 - ln_phi; // log(var - mu)
          jnll -= keep(i) * dnbinom_robust(y_i(i), s1, s2, true) * weights_i(i);
          break;
        case lognormal_family:
          jnll -= keep(i) * dlnorm(y_i(i), mu_i(i) - pow(exp(ln_phi), Type(2)) / Type(2), exp(ln_phi), true) * weights_i(i);
          break;
        case student_family:
          jnll -= keep(i) * dstudent(y_i(i), mu_i(i), exp(ln_phi), Type(3) /*df*/, true) * weights_i(i);
          break;
        case Beta_family: // Ferrari and Cribari-Neto 2004; betareg package
          s1 = mu_i(i) * exp(ln_phi);
          s2 = (Type(1) - mu_i(i)) * exp(ln_phi);
          jnll -= keep(i) * dbeta(y_i(i), s1, s2, true) * weights_i(i);
          break;
        default:
          error("Family not implemented.");
      }
    }
  }

  // ------------------ Predictions on new data --------------------------------

  if (do_predict) {
    vector<Type> proj_fe = proj_X_ij * b_j;
    // add threshold effect if specified
    if(threshold_func > 0) {
      if(threshold_func == 1) {
        // linear
        for (int i = 0; i < proj_X_ij.rows(); i++) {
          proj_fe(i) = proj_fe(i) + linear_threshold(proj_X_threshold(i), s_slope, s_cut);
        }
      } else {
        // logistic
        for (int i = 0; i < proj_X_ij.rows(); i++) {
          proj_fe(i) = proj_fe(i) + logistic_threshold(proj_X_threshold(i), s50, s95, s_max);
        }
      }
    }

    vector<Type> proj_rw_i(proj_X_ij.rows());
    for (int i = 0; i < proj_X_ij.rows(); i++) {
      proj_rw_i(i) = Type(0);
    }
    if (random_walk) {
      for (int i = 0; i < proj_X_rw_ik.rows(); i++) {
        for (int k = 0; k < proj_X_rw_ik.cols(); k++) {
          proj_rw_i(i) += proj_X_rw_ik(i, k) * b_rw_t(proj_year(i), k);
          proj_fe(i) += proj_rw_i(i);
        }
      }
    }
    REPORT(proj_rw_i);
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
    REPORT(proj_rw_i);          // random walk projections

    if (calc_se) ADREPORT(proj_eta);

    if (calc_time_totals) {
      // ------------------ Derived quantities ---------------------------------

      // Total biomass:
      vector<Type> total(n_t);
      for (int i = 0; i < proj_eta.size(); i++) {
        total(proj_year(i)) += InverseLink(proj_eta(i), link) * area_i(i);
      }
      vector<Type> log_total = log(total);
      REPORT(log_total);
      ADREPORT(log_total);

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

  if(threshold_func == 1) {
    REPORT(s_slope);
    ADREPORT(s_slope);
    REPORT(s_cut);
    ADREPORT(s_cut);
  }
  if(threshold_func == 2) {
    // report s50 and s95 for logistic function model
    REPORT(s50);
    ADREPORT(s50);
    REPORT(s95);
    ADREPORT(s95);
    REPORT(s_max);
    ADREPORT(s_max);
  }
  if (calc_quadratic_range && b_j(1) < Type(0)) {
    vector<Type> quadratic_roots = GetQuadraticRoots(b_j(1), b_j(0), Type(0.05));
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

  // jnll += nll_data + nll_omega + nll_omega_trend + nll_varphi + nll_epsilon + nll_priors;
  return jnll;
}

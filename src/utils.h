namespace sdmTMB {
template <class Type>
bool isNA(Type x) {
  return R_IsNA(asDouble(x));
}

template <class Type>
Type dcenspois(Type x, Type lambda, Type lwr, Type upr, int give_log = 0)
{
  // Should not do the obvious route due to numerical issues
  // tmp_ll = log(ppois(UPPER_i(i), mu_i(i), true) - ppois(LOWER_i(i)-1, mu_i, true));
  Type tmp_ll;
  if (lwr == upr) {  // no censorship
    tmp_ll = dpois(Type(lwr), lambda, true);
  }
  if (isNA(upr)) {  // right censored
    if (lwr == Type(0)) {
      tmp_ll = 0.0;
    }
    if (lwr > Type(0)) {
      tmp_ll = log(ppois(Type(lwr-1.0), lambda)); // F(lower-1)
      tmp_ll = logspace_sub(Type(0), tmp_ll);  // 1 - F(lower-1)
    }
  } else {
    tmp_ll = log(ppois(Type(upr), lambda)); // F(upr)
    if (lwr > Type(0)) {
      tmp_ll = logspace_sub(tmp_ll, log(ppois(Type(lwr-1.0), lambda))); // F(upr) - F(lwr-1) iff lwr>0
    }
  }
  if (give_log)
    return tmp_ll;
  else
    return exp(tmp_ll);
}

template <class Type>
Type dstudent(Type x, Type mean, Type sigma, Type df, int give_log = 0) {
  // from metRology::dt.scaled()
  // dt((x - mean)/sd, df, ncp = ncp, log = TRUE) - log(sd)
  Type logres = dt((x - mean) / sigma, df, true) - log(sigma);
  if (give_log)
    return logres;
  else
    return exp(logres);
}

template <class Type>
Type dlnorm(Type x, Type meanlog, Type sdlog, int give_log = 0) {
  Type logres = dnorm(log(x), meanlog, sdlog, true) - log(x);
  if (give_log)
    return logres;
  else
    return exp(logres);
}

// List of matrices
template <class Type>
struct LOM_t : vector<matrix<Type> > {
  LOM_t(SEXP x) {  // x = list passed from R
    (*this).resize(LENGTH(x));
    for (int i = 0; i < LENGTH(x); i++) {
      SEXP sm = VECTOR_ELT(x, i);
      (*this)(i) = asMatrix<Type>(sm);
    }
  }
};

// Function to import barrier-SPDE code
// From Olav Nikolai Breivik and Hans Skaug:
template <class Type>
struct spde_barrier_t {
  vector<Type> C0;
  vector<Type> C1;
  Eigen::SparseMatrix<Type> D0;
  Eigen::SparseMatrix<Type> D1;
  Eigen::SparseMatrix<Type> I;
  spde_barrier_t(SEXP x) { /* x = List passed from R */
    C0 = asVector<Type>(getListElement(x, "C0"));
    C1 = asVector<Type>(getListElement(x, "C1"));
    D0 = tmbutils::asSparseMatrix<Type>(getListElement(x, "D0"));
    D1 = tmbutils::asSparseMatrix<Type>(getListElement(x, "D1"));
    I = tmbutils::asSparseMatrix<Type>(getListElement(x, "I"));
  }
};

// Function to calculate Q (precision) matrix using barrier-SPDE
// From Olav Nikolai Breivik and Hans Skaug via VAST
template <class Type>
Eigen::SparseMatrix<Type> Q_spde(spde_barrier_t<Type> spde, Type kappa,
                                 vector<Type> c) {
  // using namespace Eigen;
  vector<Type> range(2);
  range(0) = sqrt(8) / kappa * c(0);
  range(1) = range(0) * c(1);

  int dimLatent = spde.D0.row(0).size();
  vector<Type> Cdiag(dimLatent);
  Eigen::SparseMatrix<Type> Cinv(dimLatent, dimLatent);

  Cdiag = spde.C0 * pow(range(0), 2) + spde.C1 * pow(range(1), 2);
  for (int i = 0; i < dimLatent; ++i) {
    Cinv.coeffRef(i, i) = 1 / Cdiag(i);
  }

  Eigen::SparseMatrix<Type> A = spde.I;
  A = A + (pow(range(0), 2) / 8) * spde.D0 + (pow(range(1), 2) / 8) * spde.D1;

  Eigen::SparseMatrix<Type> Q = A.transpose() * Cinv * A / M_PI * 2 * 3;

  return Q;
}

template <class Type>
Type minus_one_to_one(Type x) {
  return Type(2) * invlogit(x) - Type(1);
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

// FIXME no longer needed!?
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

// https://github.com/hrue/r-inla/blob/devel/r-inla.org/doc/prior/pc.matern.pdf
template <class Type>
Type pc_prior_matern(Type logtau, Type logkappa, Type matern_range,
                     Type matern_SD, Type range_prob, Type SD_prob,
                     int give_log = 0, int share_range = 0, int stan_flag = 0) {
  Type d = 2.;  // dimension
  Type dhalf = d / 2.;
  Type lam1 = -log(range_prob) * pow(matern_range, dhalf);
  Type lam2 = -log(SD_prob) / matern_SD;
  Type range = sqrt(8.) / exp(logkappa);
  Type sigma = 1. / sqrt(4. * M_PI * exp(2. * logtau) * exp(2. * logkappa));
  Type range_ll = log(dhalf) + log(lam1) + log(pow(range, -1. - dhalf)) -
                  lam1 * pow(range, -dhalf);
  Type sigma_ll = log(lam2) - lam2 * sigma;
  Type penalty = sigma_ll;
  if (!share_range) penalty += range_ll;

  // Note: these signs are + (and different from inst/jacobian-pcprior-tests)
  // because the jnll is accumulated
  if (stan_flag) {
    penalty += log(sqrt(8.)) - log(pow(range, 2.)); // P(sigma)
    Type C = sqrt(exp(lgamma(1. + dhalf)) * pow(4. * M_PI, dhalf));
    penalty += log(C) + logkappa;
  }
  // std::cout << "PC penalty: " << penalty << "\n";
  if (give_log)
    return penalty;
  else
    return exp(penalty);
}

template <class Type>
vector<Type> GetQuadraticRoots(Type a, Type b, Type threshold) {
  vector<Type> res(4);
  Type c = 1.;  // doesn't matter; setting to an arbitrary value
  Type crit_y =
      (a * pow(-b / (2. * a), 2.) + b * (-b / (2. * a)) + c) + log(threshold);
  // solve for 0 = ax2 + bx + (c - crit_y)
  c = c - crit_y;
  res(0) = -1. * (b - sqrt(pow(b, 2.) - 4. * c * a)) / (2. * a);
  res(1) = -1. * (b + sqrt(pow(b, 2.) - 4. * c * a)) / (2. * a);

  // calculate vertex of parabola
  Type xpeak = -b / (2. * a);
  // res(2) is the hi/lowpoint of parabola evaluated at xpeak
  res(2) = (a * (pow(xpeak, 2.)) + b * (xpeak) + c);

  // calculate reduction of changing from mean to +/- 1 SD
  res(3) = (a * (pow(xpeak + 1, 2.)) + b * (xpeak + 1) + c) / res(2);
  return res;
}

template <class Type>
Type linear_threshold(Type x, Type slope, Type cutpoint) {
  // linear threshold model. relationship linear up to a point then constant
  // keep all parameters unconstrained - slope and scale can be neg/pos,
  // as can cutpoint if covariate is scaled ~ N(0,1).
  Type pred;
  if (x < cutpoint) {
    pred = x * slope;
  } else {
    pred = cutpoint * slope;
  }
  return pred;
}

template <class Type>
Type logistic_threshold(Type x, Type s50, Type s95, Type scale) {
  // logistic threshold model. similar to length or size based selectivity
  // in fisheries, parameterized by the points at which f(x) = 0.5 or 0.95
  // s50 and scale are unconstrained. s95 has to be > s50 though, so modelled as
  // s95 = s50 + exp(b(1))
  // Type s95 = s50 + exp(soffset); // this done outside function
  Type pred = (scale)*Type(1.0) /
              (Type(1.0) + exp(-log(Type(19.0)) * (x - s50) / (s95 - s50)));
  return pred;
}

// from glmmTMB distrib.h
// alpha = size (dispersion param), k = truncation point, mu = mean
double rtruncated_nbinom(double alpha, int k, double mu) {
  int m;
  double mdoub;
  double p = alpha / (mu + alpha);
  double q = mu / (mu + alpha);

  if (alpha <= 0.0)
    throw std::range_error("non-positive size in k-truncated-neg-bin simulator\n");
  if (mu <= 0.0)
    throw std::range_error("non-positive mu in k-truncated-neg-bin simulator\n");
  if (k < 0)
    throw std::range_error("negative k in k-truncated-neg-bin simulator\n");

  mdoub = (k + 1.0) * p - alpha * q;
  if (mdoub < 0.0)
    mdoub = 0.0;
  m = mdoub;
  if (m < mdoub)
    m = m + 1;
  /* since p < 1.0 and q > 0.0 we have 0.0 <= mdoub < k + 1
     hence 0 <= m <= k + 1 */

  for (;;) {
    double x = rnbinom(alpha + m, p) + m;
    if (m > 0) {
      double a = 1.0;
      int j;
      double u = unif_rand();
      for (j = 0; j < m; ++j)
        a *= (k + 1 - j) / (x - j);
      if (u < a && x > k)
        return x;
    } else {
      if (x > k)
        return x;
    }
  }
}

template<class Type>
Type calc_rf_sigma(Type ln_tau, Type ln_kappa) {
  Type sigma = 1 / sqrt(Type(4.) * M_PI * exp(Type(2.) * ln_tau + Type(2.) * ln_kappa));
  return sigma;
}

// from glmmTMB distrib.h
extern "C" {
  /* See 'R-API: entry points to C-code' (Writing R-extensions) */
  double Rf_logspace_sub (double logx, double logy);
  void   Rf_pnorm_both(double x, double *cum, double *ccum, int i_tail, int log_p);
}

// from glmmTMB distrib.h
TMB_ATOMIC_VECTOR_FUNCTION(
  // ATOMIC_NAME
  logit_invcloglog
  ,
    // OUTPUT_DIM
    1,
    // ATOMIC_DOUBLE
    ty[0] = Rf_logspace_sub(exp(tx[0]), 0.);
,
  // ATOMIC_REVERSE
  px[0] = exp( logspace_add(tx[0], tx[0]-ty[0]) ) * py[0];
)

template<class Type>
Type logit_invcloglog(Type x) {
  CppAD::vector<Type> tx(1);
  tx[0] = x;
  return logit_invcloglog(tx)[0];
}

}  // namespace sdmTMB

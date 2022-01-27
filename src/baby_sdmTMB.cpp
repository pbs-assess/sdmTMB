
#include <TMB.hpp>
#include "utils.h"

using namespace density;

template<class Type>
Type objective_function<Type>::operator() () {

  using namespace density;
  // Data inputs
  DATA_VECTOR(y);
  DATA_VECTOR(pO2);
  DATA_VECTOR(temp);
  DATA_SCALAR(tref); // reference temperature for metabolic index calculation
  DATA_SCALAR(W); // assumed body weight
  DATA_INTEGER(ndata);
  DATA_MATRIX(priorsigma); // variance covariance matrix of MI parameters from bayesian fit to lab data
  DATA_VECTOR(priormeans);
  DATA_INTEGER(nprior); // posterior means of MI parameters

  // parameter inputs
  PARAMETER(Eo);
  PARAMETER(n);
  PARAMETER(logAomu);
  PARAMETER(logAosigma)
  PARAMETER(logAo);
  PARAMETER(ln_phi);
  PARAMETER(thetaf);
  PARAMETER_VECTOR(b_threshold);

  // Joint negative log-likelihood
  Type jnll = 0;

  // procedures aka transformed variables
  // Breakpoint function
  Type s_slope, s_cut;
  // these are for linear model
  s_slope = b_threshold(0);
  s_cut = b_threshold(1);
  Type kb;
  kb = 0.00008617333262145; //boltzmans constant in eV per K

  // Tweedie parameters
  Type s1, phi;
  s1 = invlogit(thetaf) + Type(1.0);
  phi = exp(ln_phi);

  // Prior parameters and calculation of Ao
  Type Ao;
  //vector<Type> parres(priormeans.size());

  Ao = exp(logAo);
  ADREPORT(Ao);
  // get prior on logdevAo
  jnll -= dnorm(logAo, logAomu, logAosigma, true);

  // calculate deviations of the vector [Eo n logAo logAosigma] from prior means

  vector<Type> parres(nprior);

  parres[0] = Eo - priormeans[0];
  parres[1] = n - priormeans[1];
  parres[2] = logAomu - priormeans[2];
  parres[3] = logAosigma - priormeans[3];

  jnll+= MVNORM(priorsigma)(parres);


  ADREPORT(phi);
  ADREPORT(s1);

  // calculate metabolic index

  vector<Type> mi(ndata);

  for(int i = 0; i<ndata; i++) {
    Type invtemp = pow(kb, -1) * (pow(temp[i], -1) - pow(tref, -1));
     mi[i] = pO2[i] * Ao * pow(W, n) * exp(Eo * invtemp);
  }


  // calculate expected value and jnll
  for(int i = 0; i<ndata; i++) {
    Type mu = exp(sdmTMB::linear_threshold(mi[i], s_slope, s_cut));
    jnll -= dtweedie(y(i), mu, phi, s1, true);
  }


  return jnll;
}


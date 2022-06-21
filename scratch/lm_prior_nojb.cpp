// Simple linear regression.
#include <TMB.hpp>
template<class Type>
  Type objective_function<Type>::operator() ()
{
  DATA_VECTOR(Y);
  DATA_VECTOR(x);
  PARAMETER(a);
  PARAMETER(b);
  PARAMETER(logSigma);
  ADREPORT(exp(2*logSigma));

  Type sigma = exp(logSigma);
  Type nll = -sum(dnorm(Y, a+b*x, exp(logSigma), true));

  nll -= dnorm(sigma, Type(0.0), Type(1), true);

  return nll;
}

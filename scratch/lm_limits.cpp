// Simple linear regression.
#include <TMB.hpp>
template<class Type>
  Type objective_function<Type>::operator() ()
{
  DATA_VECTOR(Y);
  DATA_VECTOR(x);
  PARAMETER(a);
  PARAMETER(b);
  PARAMETER(sigma);

  Type nll = -sum(dnorm(Y, a+b*x, sigma, true));
  nll -= dnorm(sigma, Type(0.0), Type(1), true);

  return nll;
}

// https://github.com/kaskr/adcomp/blob/master/tmb_examples/tweedie.cpp
// Estimating parameters in a Tweedie distribution.
#include <TMB.hpp>

template<class Type>
Type objective_function<Type>::operator() ()
{
  DATA_VECTOR(y);
  PARAMETER(mu);
  PARAMETER(phi);
  PARAMETER(p);
  Type ans = 0;
  for(int i=0; i<y.size(); i++)
    ans -= dtweedie(y(i), mu, phi, p, true);
  return ans;
}

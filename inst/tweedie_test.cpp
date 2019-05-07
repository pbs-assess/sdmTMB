// https://github.com/kaskr/adcomp/blob/master/tmb_examples/tweedie.cpp
// Estimating parameters in a Tweedie distribution.
#include <TMB.hpp>

template <class Type>
Type minus_one_to_one(Type x)
{
  return Type(2) * invlogit(x) - Type(1);
}

template<class Type>
Type objective_function<Type>::operator() ()
{
  DATA_VECTOR(y);
  PARAMETER(ln_mu);
  PARAMETER(ln_phi);
  PARAMETER(thetaf);
  Type ans = 0;
  for(int i=0; i<y.size(); i++)
    ans -= dtweedie(y(i), exp(ln_mu), exp(ln_phi), invlogit(thetaf) + Type(1.0), true);
  return ans;
}

#ifndef _NEWTON_RAPHSON_
#define _NEWTON_RAPHSON_

#include <iostream>
#include <cmath>
#include "mymacros.h"

struct NRFunct {
  // evaluate the function
  virtual real_t eval(const real_t q) const = 0;

  // evaluate the function derivative
  virtual real_t derivative(const real_t q) const =0;
};

struct NewtonRaphson
{
  constexpr auto MAXITER = (30);

  real_t solve(NRFunct& a_f, const real_t a_tol, const real_t a_guess=zero)
  {
    // initialize
    real_t s = a_guess;
    real_t err = one;
    int iter = 0;

    while (std::abs(err)>a_tol && iter<MAXITER) {
      iter++;

      const real_t f = a_f.eval(s);
      real_t df= a_f.derivative(s);
      if (std::abs(df) < tiny*f) df = SGN(df) * tiny;

      real_t ds = -f/df;
      if (ds*err<zero) ds *= half;
      s += ds;

      err = ds/std::max(s,small);
    }
    if (std::abs(err)>a_tol && iter>MAXITER)
      std::cerr << " Warning:: NewtonRapshon: solve error "
                << err << " iter " << iter << " s " << s
                << std::endl;

    return s;
  }
};

#endif

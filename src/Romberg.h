
#include "Macros.h"
#include <functional>
#include <cmath>
#include <iostream>

struct Romberg
{
  real_t integral(const real_t a, const real_t b, const real_t tol,
                  std::function<real_t (const real_t)> f, const string s="")
  {
    static const size_t ROMB_MAXITER=20;
    static const size_t ROMB_MINITER=3;
    static const size_t MAXJ=5;
    static const real_t ROMB_UNDERFLOW = 1.e-50;

    real_t h = half *(b-a);
    real_t gmax = h*(f(a)+f(b));
    real_t g[MAXJ+1];
    g[0]=gmax;

    real_t error=1.e9;
    size_t nint=1;
    size_t iter = 0;

    real_t g0;
    while ((std::abs(error)>tol || iter<ROMB_MINITER) && iter<ROMB_MAXITER)
      {
        // Calculate next trapezoidal rule approximation to integral.
        iter++;

        g0=zero;
        for (size_t k=0; k<nint; k++) {
          g0 += f(a+(two*k+1)*h);
        }

        g0 = half*g[0]+h*g0;
        h *= half;
        nint *= 2;
        size_t jmax = iter<MAXJ ? iter : MAXJ;
        real_t fourj=one;

        for (size_t j=0; j<jmax; j++) {
          // Use Richardson extrapolation.
          fourj=four*fourj;
          real_t g1=g0+(g0-g[j])/(fourj-one);
          g[j]=g0;
          g0=g1;
        }

        error = std::abs(g0)>ROMB_UNDERFLOW ? one-gmax/g0 : gmax;

        gmax=g0;
        g[jmax]=g0;
        //    if (i>5 && std::abs(error)< tol) return g0;
      }
    if (error>tol) {
      std::cerr << "Rombint on " << s << " I="<<g0
                << ", a="<<a<<", b="<<b
                << " ... failed to converge: error : " << error << " ! tol="<< tol <<'\n';
      //exit(0);
    }

    return g0;
  }
};


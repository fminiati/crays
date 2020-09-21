//
// Copyright (C) 2020 Francesco Miniati <francesco.miniati@gmail.com>
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in all
// copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
// SOFTWARE.
//
#ifndef ROMBERG_H
#define ROMBERG_H

#include "Macros.h"
#include <functional>
#include <cmath>
#include <iostream>

// perform Romberg's integration;
// implementation based on Numerical Recipes' Fortran77 version.
namespace fm {
    struct Romberg
    {
        real_t integral(const real_t a, const real_t b, const real_t tol,
                        std::function<real_t(const real_t)> f, const std::string &&s = {})
        {
            static const size_t ROMB_MAXITER = 20;
            static const size_t ROMB_MINITER = 3;
            static const size_t MAXJ = 5;
            static const real_t ROMB_UNDERFLOW = 1.e-50;

            real_t h = half * (b - a);
            real_t gmax = h * (f(a) + f(b));
            real_t g[MAXJ + 1];
            g[0] = gmax;

            real_t error = 1.e9;
            size_t nint = 1;
            size_t iter = 0;

            real_t g0;
            while ((std::abs(error) > tol || iter < ROMB_MINITER) && iter < ROMB_MAXITER)
            {
                // Calculate next trapezoidal rule approximation to integral.
                iter++;

                g0 = zero;
                for (size_t k = 0; k < nint; k++)
                {
                    g0 += f(a + (two * k + 1) * h);
                }

                g0 = half * g[0] + h * g0;
                h *= half;
                nint *= 2;
                size_t jmax = iter < MAXJ ? iter : MAXJ;
                real_t fourj = one;

                for (size_t j = 0; j < jmax; j++)
                {
                    // apply Richardson's extrapolation.
                    fourj = four * fourj;
                    real_t g1 = g0 + (g0 - g[j]) / (fourj - one);
                    g[j] = g0;
                    g0 = g1;
                }

                error = std::abs(g0) > ROMB_UNDERFLOW ? one - gmax / g0 : gmax;

                gmax = g0;
                g[jmax] = g0;
            }
            if (error > tol)
            {
                std::cerr << "Rombint on " << s << " I=" << g0
                          << ", a=" << a << ", b=" << b
                          << " ... failed to converge: error : " << error << " ! tol=" << tol << '\n';
                //exit(0);
            }

            return g0;
        }
    };
}; // namespace fm
#endif

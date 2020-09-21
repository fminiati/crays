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
#ifndef NEWTONRAPHSON_H
#define NEWTONRAPHSON_H

#include <cmath>
#include <iostream>

namespace fm {
    struct NRFunct
    {
        // evaluate the function
        virtual double eval(const double q) const = 0;

        // evaluate the function derivative
        virtual double derivative(const double q) const = 0;
    };

    struct NewtonRaphson
    {
        static constexpr double zero = 0.0;
        static constexpr double half = 0.5;
        static constexpr double one = 1.0;
        static constexpr double small = 1.e-6;
        static constexpr double tiny = 1.e-9;
        static constexpr int MAXITER = 30;

        double solve(NRFunct &a_f, const double a_tol, const double a_guess = zero)
        {
            // initialize
            double s = a_guess;
            double err = one;
            int iter = 0;

            while (std::abs(err) > a_tol && iter < MAXITER)
            {
                iter++;

                const double f = a_f.eval(s);
                double df = a_f.derivative(s);
                if (std::abs(df) < tiny * f)
                    df = (df >= 0 ? 1 : -1) * tiny;

                double ds = -f / df;
                if (ds * err < zero)
                    ds *= half;
                s += ds;

                err = ds / std::max(s, small);
            }
            if (std::abs(err) > a_tol && iter > MAXITER)
                std::cerr << " Warning:: NewtonRapshon: solve error "
                          << err << " iter " << iter << " s " << s
                          << std::endl;

            return s;
        }
    };
}; // namespace fm
#endif

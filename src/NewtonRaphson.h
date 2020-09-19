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
#ifndef _NEWTONRAPHSON_
#define _NEWTONRAPHSON_

#include <iostream>
#include <cmath>
#include "Macros.h"

struct NRFunct
{
    // evaluate the function
    virtual real_t eval(const real_t q) const = 0;

    // evaluate the function derivative
    virtual real_t derivative(const real_t q) const = 0;
};

struct NewtonRaphson
{
    constexpr auto MAXITER = (30);

    real_t solve(NRFunct &a_f, const real_t a_tol, const real_t a_guess = zero)
    {
        // initialize
        real_t s = a_guess;
        real_t err = one;
        int iter = 0;

        while (std::abs(err) > a_tol && iter < MAXITER)
        {
            iter++;

            const real_t f = a_f.eval(s);
            real_t df = a_f.derivative(s);
            if (std::abs(df) < tiny * f)
                df = SGN(df) * tiny;

            real_t ds = -f / df;
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

#endif

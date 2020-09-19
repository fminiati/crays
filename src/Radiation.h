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
#ifndef _RADIATION_
#define _RADIATION_

#include <cmath>
#include <chrono>
#include <functional>
#include <string>
#include "Macros.h"
#include "Romberg.h"
#include "PhysicalConstants.h"

namespace Radiation
{
    static constexpr real_t _TOL = 5.e-3;

    struct EmissionProc
    {
        EmissionProc() = default;
        EmissionProc(EmissionProc const &) = default;
        EmissionProc(const real_t pmin, const real_t pmax,
                     std::function<real_t(const real_t)> f)
            : _pmin(pmin), _pmax(pmax), _fp(f)
        {
        }
        real_t _pmin;
        real_t _pmax;
        std::function<real_t(const real_t)> _fp;
    };

    struct Synchrotron : public EmissionProc
    {
        Synchrotron() = default;
        Synchrotron(Synchrotron const &) = default;
        Synchrotron(const real_t pmin, const real_t pmax,
                    std::function<real_t(const real_t)> f)
            : EmissionProc(pmin, pmax, f)
        {
        }

        // compute synchrotron spectrum at nu scaled by critical frequency
        // nc=1.5(eB/mc)
        void operator()(std::vector<real_t> &a_spectrum,
                        const std::vector<real_t> &a_nu, const size_t a_beg,
                        const size_t a_stride, const size_t a_size) const
        {
            for (size_t i = a_beg; i < a_size; i += a_stride)
            {
                const real_t nu = a_nu[i];
                // integrand function in log space
                auto ksync = [&_f = _fp, _nu = nu, _plo = _pmin, _phi = _pmax](const real_t a_lx)
                {
                    Romberg r;
                    const real_t x = exp(a_lx);

                    // factor due to pitch angle average
                    const real_t mulo = sqrt(one - std::min(one, pow(_nu / (x * pow(_plo, 2)), 2)));
                    const real_t muhi = sqrt(one - std::min(one, pow(_nu / (x * pow(_phi, 2)), 2)));
                    if (mulo > muhi)
                        return zero;

                    auto fsina = [&_f = _f, _nx = _nu / x](const real_t mu) {
                        // from x= (nu/nc) p^-2; get p[mec]
                        const real_t sina = std::max(tiny, sqrt(one - mu * mu));
                        const real_t p = sqrt(_nx / sina);
                        return (_f(p) / sqrt(sina));
                    };
                    const real_t muaf = r.integral(zero, one, 10 * _TOL, fsina, "fsina");
                    //const real_t muaf= r.integral(mulo,muhi,10*_TOL,fsina,"fsina");

                    // M.Kh.Khokonov. JETP, V.99, No.4, pp. 690-707 (2004)
                    auto k53 = [_x = x](const real_t lz) {
                        static const real_t sqrt3 = sqrt(three);
                        const real_t z = exp(lz);
                        const real_t z23 = third * z * z;
                        const real_t b = (1 + 4 * z23) * sqrt(1 + z23);
                        const real_t a = (1 + 4 * z23 * (3 + 4 * z23)) / b;
                        // multiply by z to compensate dlz
                        return (z * sqrt3 * a / b * exp(-b * _x));
                    };
                    const real_t fsy = x * r.integral(-hundred, three, _TOL, k53, "k53");
                    const real_t x5h = exp(-2.5 * a_lx);

                    // altern. 1.8*x^1/3*e^-x (AM Thompson A&A 1990, 240, 209-215)
                    //static constexpr real_t b= 7.e0/8e0;
                    //static constexpr real_t c=11.e0/8e0;
                    //const real_t fsy= 0.87e0 *exp((third-2.5)*a_lx -c*exp(b*a_lx));

                    // multiply by x to compensate dlx
                    return (x * fsy * x5h * muaf);
                };

                // integrate in log scale
                Romberg r;
                a_spectrum[i] = pow(nu, 1.5) * r.integral(-hundred, three, _TOL, ksync, "ksync");
            }
        }
    };
}; // namespace Radiation

#endif

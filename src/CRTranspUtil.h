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
#ifndef CRTRANSPUTIL_H
#define CRTRANSPUTIL_H

#include <vector>
#include <iomanip>
#include <cassert>
#include <map>
#include <cmath>
#include "Macros.h"

namespace cr_transport
{
    // time dependent model of plasma:
    // n=number density, T=temperature, B=magnetic field
    // divv=average divergence of turbulent velocity on scale _ell
    // zeta=power-law index of dv scaling with separation: dv(l)~dv(l/_ell)^z
    // t=time
    struct PlasmaModel
    {
        PlasmaModel()
            : _it(0) {}

        ~PlasmaModel() {}

        // data
        real_t _ell;
        std::vector<real_t> _n, _T, _B, _divv, _zeta, _t;

        // return variable at time t
        real_t get(const real_t a_t, const std::string a_var)
        {
            assert(a_t >= _t.front() && a_t <= _t.back());

            update_index(a_t);

            auto t_interp = [i = _it, t = _t](const auto ta, const auto &f) {
                return ((f[i] - f[i - 1]) / (t[i] - t[i - 1]) * (ta - t[i]) + f[i]);
            };

            if (a_var == "density")
                return t_interp(a_t, _n);
            else if (a_var == "temperature")
                return t_interp(a_t, _T);
            else if (a_var == "mag-field")
                return t_interp(a_t, _B);
            else if (a_var == "divv-turb")
                return t_interp(a_t, _divv);
            else if (a_var == "zeta-turb")
                return t_interp(a_t, _zeta);

            return zero;
        }

    private:
        // iterpolation index
        void update_index(const real_t a_t)
        {
            assert(_t.size() > 0);

            if (!(_it > 0 && a_t > _t[_it - 1] && a_t < _t[_it]))
                while (a_t >= _t[_it])
                    ++_it;
        }

    private:
        size_t _it;
    };

    // node for lagrangian description of phase space distribution
    // function
    struct DistFunctNode
    {
        DistFunctNode() = default;
        DistFunctNode(DistFunctNode const &) = default;
        //DistFunctNode(DistFunctNode &&)=default;

        real_t f;
        real_t p;
    };

    struct Injection
    {
        Injection() = default;
        Injection(Injection const &) = default;
        virtual ~Injection() {}

        // actual function
        virtual real_t operator()(const real_t x) const = 0;
        // low and high momentum cut-off
        real_t _plo, _phi;
    };

    // f(p)=1/(pi*Gamma) / (1+(p-p0)^2/Gamma^n), normalise to
    // injection time
    struct GammaInjection : public Injection
    {
        GammaInjection() = default;
        GammaInjection(GammaInjection const &) = default;
        virtual ~GammaInjection() {}

        // actual function
        virtual real_t operator()(const real_t p) const
        {
            return one / (Pi * _Gamma * _dt) / (one + std::pow((p - _p0) / _Gamma, _idx));
        }

        size_t _idx;
        real_t _p0, _Gamma, _dt;
    };

    template <typename I>
    struct InjectionOp
    {
        InjectionOp() = default;
        InjectionOp(const real_t a_pmin,
                    const real_t a_dlp,
                    const I &a_i)
            : _pmin(a_pmin), _dlp(a_dlp), _inject(a_i)
        {
        }

        //
        void operator()(std::vector<DistFunctNode> &a_j,
                        const real_t a_dt, const size_t a_beg,
                        const size_t a_stride, const size_t a_size)
        {
            for (size_t b = a_beg; b < a_size; b += a_stride)
            {
                const real_t p = _pmin * std::exp(b * _dlp);
                a_j[b].p = p;
                if (p >= _inject._plo && p <= _inject._phi)
                    a_j[b].f = a_dt * _inject(p);
            }
        }

    protected:
        // p-grid quantities
        real_t _pmin, _dlp;
        // injection funtion: returns f(p)
        I _inject;
    };
}; // namespace cr_transport

#endif

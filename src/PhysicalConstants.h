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
#ifndef _PHYS_CONSTANTS_
#define _PHYS_CONSTANTS_

#include "Macros.H"

namespace PhysConstCgs
{
    static constexpr real_t c = 2.99792458e10;    // speed of light
    static constexpr real_t h_p = 6.62606936e-27; // Planck constant
    static constexpr real_t e = 4.80320441e-10;   // electron charge
    static constexpr real_t m_e = 9.10938262e-28; // electron mass
    static constexpr real_t m_p = 1.67262171e-24; // proton mass

    static constexpr real_t G_N = 6.6742e-8;           // Newton constant
    static constexpr real_t G_F = 1.6637e-5;           // Fermi  constant
    static constexpr real_t k_B = 1.3806505e-16;       // Boltzmann constant
    static constexpr real_t sigma_T = 6.652458558e-25; // Thomson cross section

    static constexpr real_t rho_c = 1.87837e-29;    // rho critical/h^2
    static constexpr real_t H_0 = 3.2407797216e-18; // Hubble const/h
    static constexpr real_t M_sun = 1.98844e33;     // Sun mass
    static constexpr real_t L_sun = 3.846e33;       // Sun luminosity

    static constexpr real_t AU = 1.49597870660e13; // Astronomical Unit
    static constexpr real_t yr = 3.15569252e7;     // Year
    static constexpr real_t km = 1.e5;             // km
    static constexpr real_t pc = 3.0856775807e18;  // Parsec
    static constexpr real_t kpc = 3.0856775807e21; // KiloParsec
    static constexpr real_t Mpc = 3.0856775807e24; // MegaParsec

    static constexpr real_t T_CMB = 2.725e0;     // CMB Temperature  @z=0
    static constexpr real_t U_CMB = 4.17175e-13; // CMB Energy dens. @z=0
    static constexpr real_t n_CMB = 4.1050e2;    // CMB ph. # dens. @z=0
};                                               // namespace PhysConstCgs

#endif

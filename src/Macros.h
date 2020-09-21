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
#ifndef MACROS_H
#define MACROS_H

// my variables
using real_t = double;

template <typename T>
int SGN(const T t) {return (t>T(0) ? 1 : -1);};
namespace fm {
    static constexpr auto zero = (0.0e0);
    static constexpr auto half = (0.5e0);
    static constexpr auto one = (1.0e0);
    static constexpr auto two = (2.0e0);
    static constexpr auto three = (3.0e0);
    static constexpr auto four = (4.0e0);
    static constexpr auto five = (5.0e0);
    static constexpr auto six = (6.0e0);
    static constexpr auto seven = (7.0e0);
    static constexpr auto eight = (8.0e0);
    static constexpr auto nine = (9.0e0);
    static constexpr auto ten = (10.0e0);
    static constexpr auto twenty = (20.0e0);
    static constexpr auto tenth = (0.100e0);
    static constexpr auto eighth = (0.125e0);
    static constexpr auto fifth = (0.200e0);
    static constexpr auto fourth = (0.250e0);
    static constexpr auto third = (one / three);
    static constexpr auto Pi = (3.14159265358979323846e0);
    static constexpr auto small = (1.e-6);
    static constexpr auto tiny = (1.e-9);
    static constexpr auto hundred = (1.0e2);
    static constexpr auto huge = (1.0e100);
}; // namespace fm

#endif

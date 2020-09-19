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
#include <cmath>
#include <iostream>
#include <sstream>
#include <chrono>
#include "FileParser.h"
#include "PhysicalConstants.h"
#include "CRTransport.h"
#include "CRTranspUtil.h"
#include "PThreadUtil.h"
#include "Radiation.h"

using namespace Radiation;

template <typename C>
void plot_histo(const C &c, const size_t n, const std::string&& s = {})
{
    std::cout << s << '\n';
    real_t maxv{};
    for (auto v : c)
        maxv = std::max(maxv, v.f);
    for (auto v : c)
    {
        const int x = static_cast<int>(n * v.f / maxv);
        if (x >= 1)
            std::cout << std::fixed
                      << std::setfill('-')
                      << std::setprecision(2)
                      << v.p << " " << std::string(x, '*') << '\n';
    }
};


namespace cr_transport
{
    using dfnode_t = std::vector<DistFunctNode>;

    void CRTransport::setup(const std::string a_file)
    {
        //
        std::cout << " setup begins " << '\n';

        // read input file
        FileParser input(a_file);

        // number of threds
        input.get_item(_num_threads, "crt.num_threads");

        real_t unitTime = one;
        input.get_item(unitTime, "crt.unit_of_time[Gyr]");
        _unit_time_sec = unitTime * 1.e9 * PhysConstCgs::yr;

        // time
        input.get_item(_t_init, "crt.initial_time[Gyr]");
        input.get_item(_t_end, "crt.final_time[Gyr]");
        // normalise
        _t_init /= unitTime;
        _t_end /= unitTime;
        input.get_item(_max_step, "crt.max_num_steps");

        // cfl number for time integration
        input.get_item(_cfl, "crt.CFL_number");
        assert(_cfl > zero && _cfl < one);

        // output
        input.get_item(_output_filename, "crt.output_file");
        input.get_item(_output_sync_spectrum, "crt.output_sync_spectrum");
        input.get_item(_dt_output, "crt.t_output_interval[Gyr]");
        // initialise t-output
        _t_output = zero;

        {   // compute synchtrotron emission around critical sync frequency
            input.get_item(_num_nu, "crt.sync.num_nusync");
            input.get_item(_numin, "crt.sync.nu_min[nu_crit]");
            input.get_item(_numax, "crt.sync.nu_max[nu_crit]");
            _dlnu = log(_numax / _numin) / (_num_nu - 1);
            _nu.resize(_num_nu);
            for (size_t n = 0; n < _num_nu; ++n)
                _nu[n] = _numin * exp(n * _dlnu);
        }

        {   // icm model
            // dens must be in cm^-3 and B in G
            int n_steps = 0;
            input.get_item(n_steps, "crt.icm_model.num_steps");
            input.get_items(_icm._t, n_steps, "crt.icm_model.time[sec]");
            input.get_items(_icm._n, n_steps, "crt.icm_model.gas_num_dens[cm^-3]");
            input.get_items(_icm._T, n_steps, "crt.icm_model.gas_temperature[keV]");
            input.get_items(_icm._B, n_steps, "crt.icm_model.magnetic_field[G]");
            input.get_items(_icm._divv, n_steps, "crt.icm_model.dv_turb_at_ell[km/sec]");
            input.get_item(_icm._ell, "crt.icm_model.ell_of_dv_turb[kpc]");
            input.get_items(_icm._zeta, n_steps, "crt.icm_model.pwrlaw_idx_dv_turb");

            // normalise
            _icm._ell *= PhysConstCgs::kpc;
            for (auto &t : _icm._t)
                t /= _unit_time_sec;
            std::cout << "times: ";
            for (auto &t : _icm._t)
                std::cout << t << "  ";
        }

        {   // transport parameters
            input.get_item(_lmfp, "crt.cr_model.mfp_at_p0[kpc]");
            input.get_item(_pmfp, "crt.cr_model.p0_of_mfp[mc]");
            input.get_item(_eta, "crt.cr_model.pwrlaw_idx_mfp");

            // normalise
            _lmfp *= PhysConstCgs::kpc;
            // flight time across mean-free-path
            _tau_mfp = _lmfp / PhysConstCgs::c / _unit_time_sec;
        }

        {   // rescale turb velocity increment to mfp scale and convert to <div(v)>
            // normalization factor
            const real_t A = PhysConstCgs::km / _lmfp * _unit_time_sec;
            const real_t x = log(_lmfp / _icm._ell);
            for (size_t t = 0; t < _icm._divv.size(); ++t)
            {
                _icm._divv[t] *= A * std::tgamma(_icm._zeta[t]) * exp(x * _icm._zeta[t]);
            }
        }

        {   // grid parameters
            input.get_item(_num_nodes, "crt.fc.num_nodes");
            input.get_item(_pmin, "crt.fc.p_minimum[mc]");
            input.get_item(_pmax, "crt.fc.p_maximum[mc]");
        }

        {   // inject pseudo-lorentzian with input momentum-width within time-window
            input.get_item(_inj._p0, "crt.inj.p0_injection[mc]");
            input.get_item(_inj._Gamma, "crt.inj.Gamma_injection[mc]");
            input.get_item(_inj._idx, "crt.inj.index_injection");

            // time-window
            input.get_item(_ti_inj, "crt.inj.time_injection_start[Gyr]");
            input.get_item(_te_inj, "crt.inj.time_injection_end[Gyr]");
            assert(_te_inj > _ti_inj);
            // normalise
            _ti_inj /= unitTime;
            _te_inj /= unitTime;

            // injection rate
            _inj._dt = (_te_inj - _ti_inj);

            input.get_item(_inj_multiplicity, "crt.inj.inj_multiplicity");
            // max range of fc...
            real_t max_range = one;
            input.get_item(max_range, "crt.inj.max_fc_range");
            _inj._plo = std::max(_pmin, _inj._p0 - pow(_inj._Gamma * max_range, one / _inj._idx));
            _inj._phi = std::min(_pmax, _inj._p0 + pow(_inj._Gamma * max_range, one / _inj._idx));
        }

        // coefficients for lepton loss rates in s^-1 (Coulomb+bremsstrahlung+Compton):
        // -dlnp/dt= aa+aa2*lg(n) + bb*p + cc*p^2

        {   // set timestep
            const real_t nicm = _icm.get(_t_init, "density");
            const real_t divv = _icm.get(_t_init, "divv-turb");
            const real_t zeta = _icm.get(_t_init, "zeta-turb");
            //coeff not sensitive to p, use geometric mean
            const real_t plo = _inj._plo;
            const real_t phi = _inj._phi;
            const real_t pav = sqrt(plo * phi);

            const real_t gmma = sqrt(one + pav * pav);
            const real_t beta = pav / gmma;
            const real_t aa = 1.50e-14 * (73.56e0 + log(gmma)) / (beta * beta) * _unit_time_sec;
            const real_t bb = 1.395e-16 * (log(two * gmma) - one / three) / (beta * beta) * _unit_time_sec;
            const real_t cc = 3.28e-8 * PhysConstCgs::U_CMB / beta * _unit_time_sec;
            const real_t aa2 = -1.50e-14 * _unit_time_sec;

            // transport coeff
            const real_t Dp = divv * exp(_eta * (zeta - one) * log(sqrt(plo * phi) / _pmfp));
            const real_t ae = (aa + aa2 * log(nicm)) * nicm;
            const real_t be = bb * nicm + Dp;
            const real_t clng_rate = std::max(abs(ae / plo + be + cc * plo), abs(ae / phi + be + cc * phi));

            // initial timestep is a small fraction of min cooling time
            const real_t dlp = log(phi / plo) / _num_nodes;
            _new_dt = 0.1 * _cfl / clng_rate * (exp(dlp) - one);
            assert(_new_dt > zero);

            std::cout << " paramenters:" << '\n';
            std::cout << "        t_init:" << _t_init << '\n';
            std::cout << "         t_end:" << _t_end << '\n';
            std::cout << "     momentum-grid:" << '\n';
            std::cout << "       n-nodes:" << _num_nodes << '\n';
            std::cout << "     transport-coeff:" << '\n';
            std::cout << "         l_mfp:" << _lmfp << '\n';
            std::cout << "         p_mfp:" << _pmfp << '\n';
            std::cout << "           eta:" << _eta << '\n';
            std::cout << "           tau:" << _tau_mfp << '\n';
            std::cout << "         injection:" << '\n';
            std::cout << "           p0:" << _inj._p0 << '\n';
            std::cout << "        Gamma:" << _inj._Gamma << '\n';
            std::cout << "        index:" << _inj._idx << '\n';
            std::cout << "       p_low :" << _inj._plo << '\n';
            std::cout << "       p-high:" << _inj._phi << '\n';
            std::cout << "       ti_inj:" << _ti_inj << '\n';
            std::cout << "       te_inj:" << _te_inj << '\n';
            std::cout << "             other:" << '\n';
            std::cout << "          clf:" << _cfl << '\n';
            std::cout << " " << '\n';
            std::cout << " f-injection" << '\n';

            // for (size_t ip=0; ip<_num_nodes; ip++) {
            //   const real_t p=plo * exp(ip*dlp);
            //   const int x=static_cast<int>(50*_inj(p)/_inj(plo));
            //   if (x>=1) std::cout << std::fixed << p << " " << std::string(x,'*') << '\n';
            // }
        }

        // write output
        output(zero, 0);

        std::cout << " setup ends " << '\n';
    }

    // lagrangian scheme: advect f in momentum space
    void CRTranspOp::operator()(std::vector<DistFunctNode> &a_fc,
                                real_t &a_dt,
                                const size_t a_beg,
                                const size_t a_std,
                                const size_t a_size) const
    {
        // bin and node centered slopes: qb(p), qn(p)
        std::vector<real_t> qn(a_size, zero);
        {
            size_t ip = a_beg;
            if (ip == 0)
            {
                qn[ip] = -log(a_fc[ip + 1].f / a_fc[ip].f) / log(a_fc[ip + 1].p / a_fc[ip].p);
                ip += a_std;
            }
            while (ip < a_fc.size() - 1)
            {
                qn[ip] = -log(a_fc[ip + 1].f / a_fc[ip - 1].f) / log(a_fc[ip + 1].p / a_fc[ip - 1].p);
                ip += a_std;
            }
            if (ip == a_fc.size() - 1)
            {
                qn[ip] = -log(a_fc[ip].f / a_fc[ip - 1].f) / log(a_fc[ip].p / a_fc[ip - 1].p);
            }
        }

        {
            // min_dt is used to set next timestep
            real_t dt_min = huge;

            // total energy density with B in units of CMB energy density
            const real_t Bcmb = sqrt(8 * Pi * PhysConstCgs::U_CMB);
            const real_t Ue = one + _B * _B / (Bcmb * Bcmb);
            const real_t logn = log(_n);

            for (size_t ip = a_beg; ip < a_size; ip += a_std)
            {
                // current momentum
                const real_t p = a_fc[ip].p;

                // CMB energy density at z=0
                constexpr auto UCMB = PhysConstCgs::U_CMB;

                const real_t gmma = sqrt(one + p * p);
                const real_t beta = p / gmma;
                const real_t aa = 1.50e-14 * (73.56e0 + log(gmma)) / (beta * beta) * _utime;
                const real_t bb = 1.395e-16 * (log(two * gmma) - one / three) / (beta * beta) * _utime;
                const real_t cc = 3.28e-8 * UCMB / beta * _utime;
                const real_t aa2 = -1.50e-14 * _utime;

                // transport coefficients: Dp=<divv>, hDpp=Dpp/p^2
                const real_t Dp = _divv * exp(_eta * (_zeta - one) * log(p / _pmfp));
                const real_t hDpp = Dp * Dp * exp(_eta * log(p / _pmfp)) * _tau_mfp;

                // defines the total energy loss/gain coefficients b=a+bp+cp**2
                // Dp is Fermi-I and q*Dpp is Fermi-II acceleration respectively
                const real_t ae = (aa + aa2 * logn) * _n;
                const real_t be = bb * _n + Dp - qn[ip] * hDpp;
                const real_t ce = cc * Ue;
                // find pnew defined by: dt=\int_p0^pt dp/pdot
                real_t z;
                {
                    const real_t delta = be * be - four * ae * ce;
                    const real_t dbdp = be + two * ce * p;
                    if (delta < -tiny)
                    {
                        const real_t sqrd = sqrt(-delta);
                        z = atan(dbdp / sqrd) - half * sqrd * a_dt;
                        z = tan(z) * sqrd;
                    }
                    else if (delta > tiny)
                    {
                        const real_t sqrd = sqrt(delta);
                        z = exp(-sqrd * a_dt) * (sqrd - dbdp) / (sqrd + dbdp);
                        z = sqrd * (one - z) / (one + z);
                    }
                    else
                    {
                        z = dbdp / (one + half * dbdp * a_dt);
                    }
                }
                const real_t pnew = (half * (z - be) / ce);

                if (pnew <= zero || pnew / pnew != one)
                    std::cout << "ip=" << ip << ", p=" << p << ", q=" << qn[ip] << ", be=" << be << ", Dp=" << Dp
                              << " hDpp=" << hDpp << ", n=" << _n << ", bb=" << bb
                              << " beta=" << beta << ", g=" << gmma << ", pnew=" << pnew << '\n';
                assert(pnew > zero);

                // p derivative of pdot
                const real_t dpdotdp = -(be + ce * (p + pnew));
                // update p and fp
                a_fc[ip].f *= pow(p / pnew, two) * exp(-a_dt * dpdotdp);
                a_fc[ip].p = pnew;

                // update max cooling rate
                const real_t pdot = -(ae + half * (p + pnew) * (be + ce * half * (p + pnew)));
                const long long int ips = ip + SGN(pdot);
                if (ips >= 0 && ips < a_fc.size())
                {
                    const real_t dp = _cfl * (a_fc[ips].p - a_fc[ip].p);
                    // if (abs(dp/(pdot+dp*dpdotdp))<dt_min) {
                    //   if (dt_min<5.e-6)
                    //   //if (a_beg==0)
                    //   std::cout << "ip="<<ip<<", ips="<<ips<<", dtmin="<<dt_min
                    //             <<", p="<<pnew<<", pdot="<<pdot<<", dp="<<dp
                    //             <<", dp*dpdot="<< dp*dpdotdp<<'\n';
                    // }
                    dt_min = std::min(dt_min, abs(dp / (pdot + dp * dpdotdp)));
                }
            }
            // use a_dt to return dt_min
            a_dt = dt_min;
        }
    }

    // conservative scheme: flux update
    void CRTransport::injection(std::vector<DistFunctNode> &a_f, const real_t a_dt)
    {
        std::cout << "    injection(beg): ";

        // rebin after injection
        // set min and max...
        const real_t pmin = (!a_f.empty() ? std::min(a_f.front().p, _inj._plo) : _inj._plo);
        const real_t pmax = (!a_f.empty() ? std::max(a_f.back().p, _inj._phi) : _inj._phi);
        const real_t dlp = std::log(pmax / pmin) / (_num_nodes - 1);

        // injection operator
        InjectionOp<GammaInjection> injOp(pmin, dlp, _inj);
        std::vector<dfnode_t> sinj(_num_threads, dfnode_t(_num_nodes));
        ThreadUt<InjectionOp<GammaInjection>, dfnode_t, const real_t> th_ut{_num_threads,
                                                                            _num_nodes};
        th_ut.launch_threads(injOp, sinj, a_dt);

        {
            auto t_beg = std::chrono::high_resolution_clock::now();
            // new distribution
            std::vector<DistFunctNode> fnew(_num_nodes, DistFunctNode{zero, zero});
            // add injection
            for (size_t t = 0; t < _num_threads; ++t)
            {
                for (size_t i = t; i < _num_nodes; i += _num_threads)
                {
                    fnew[i].p = sinj[t][i].p;
                    fnew[i].f = sinj[t][i].f;
                }
            }
            // add f by remapping
            constexpr real_t eps = 5.e-13;
            for (const auto &f : a_f)
            {
                const size_t b = static_cast<size_t>(std::floor(log(f.p / pmin) / dlp));
                if (b + 1 < _num_nodes)
                {
                    const real_t w = (f.p - fnew[b].p) / (fnew[b + 1].p - fnew[b].p);
                    //std::cout <<"b="<<b<<std::setprecision(18)<<", w="<<w<<", nodes="<<_num_nodes<<'\n';
                    assert(w >= zero - eps && w <= one + eps);
                    fnew[b].f += (one - w) * f.f;
                    fnew[b + 1].f += w * f.f;
                }
                else if (b + 1 == _num_nodes)
                    fnew[b].f += f.f;
                else
                {
                    std::cout << " error: b=" << b << '\n';
                    exit(0);
                }
            }
            {   // remove empty elements
                size_t nsz = fnew.size();
                for (auto n = fnew.begin(); n != fnew.end();)
                {
                    if (n->f == zero)
                        n = fnew.erase(n);
                    else
                        ++n;
                }
                nsz -= fnew.size();
                if (nsz > 0)
                    std::cout << " removed " << nsz << " particles" << '\n';
            }
            // swap old for new
            a_f.swap(fnew);

            auto t_end = std::chrono::high_resolution_clock::now();
            _timings["merging"] += t_end - t_beg;
        }

        std::cout << "    injection(end): ";
    }

    // conservative scheme: flux update
    size_t CRTransport::advance(const real_t a_t, const real_t a_dt)
    {
        std::cout << "  advance(beg): " << std::scientific << "by dt=" << a_dt << ", at t=" << a_t
                  << '\n';

        // temp distr function
        std::vector<DistFunctNode> &f = _fc;

        // source predictor step: injection window
        if (a_t >= _ti_inj && a_t <= _te_inj)
        {
            auto t_beg = std::chrono::high_resolution_clock::now();
            injection(f, half * a_dt);
            auto t_end = std::chrono::high_resolution_clock::now();
            _timings["injection"] += t_end - t_beg;
        }

        // dt_min sets the next timestep
        real_t dt_min = huge;

        {   // thread this time consuming spectrum calculation
            // create pusher
            CRTranspOp crtOp(_icm.get(a_t, "density"),
                             _icm.get(a_t, "temperature"),
                             _icm.get(a_t, "mag-field"),
                             _icm.get(a_t, "divv-turb"),
                             _icm.get(a_t, "zeta-turb"),
                             _lmfp, _pmfp, _eta, _tau_mfp,
                             _unit_time_sec,
                             _cfl);
            // args to be passed to sync function
            std::vector<real_t> dts(_num_threads, a_dt);

            auto t_beg = std::chrono::high_resolution_clock::now();
            ThreadUt<const CRTranspOp, dfnode_t, real_t> th_ut{_num_threads, f.size()};
            th_ut.launch_threads(crtOp, f, dts);
            auto t_end = std::chrono::high_resolution_clock::now();
            _timings["cr_transp"] += t_end - t_beg;

            // get min dt
            for (const auto dt : dts)
                dt_min = std::min(dt_min, dt);
        }

        // injection window, we assume the non-stiff regime
        if (a_t >= _ti_inj && a_t <= _te_inj)
        {
            auto t_beg = std::chrono::high_resolution_clock::now();
            injection(f, half * a_dt);
            auto t_end = std::chrono::high_resolution_clock::now();
            _timings["injection"] += t_end - t_beg;
        }

        {   // boundary conditions
            auto t_beg = std::chrono::high_resolution_clock::now();
            for (auto x = f.begin(); x != f.end();)
            {
                if (x->p >= _pmin && x->p <= _pmax)
                    ++x;
                else
                    x = f.erase(x);
            }
            auto t_end = std::chrono::high_resolution_clock::now();
            _timings["bc"] += t_end - t_beg;
        }
        // copy back: fc <- f
        //  _fc=f;

        // new timestep
        _new_dt = dt_min;

        std::cout << "  advance(end): " << '\n';

        return f.size();
    }

    //
    void CRTransport::output(const real_t a_t, const unsigned a_step)
    {
        if (a_step == 0)
        {
            std::cout << "  output(beg): " << '\n';

            {   // open file, pring header, momentum grid and initial function
                std::ofstream file(_output_filename, std::ios_base::trunc);
                file << "# time evolution of distribution function" << '\n';
                file << "# first  line: num of nodes" << '\n';
                file << "# second line: initial time, step, p-nodes, f values at p-nodes" << '\n';
                file << "# ...2+nth line: nth time, step, p-nodes  and corresponding f(p)" << '\n';
                file << "# injection is active from: " << _ti_inj << " --to-- " << _te_inj << " Gyr" << '\n';
                file << "# nbins=" << _num_nodes << '\n';
                file.close();
            }

            if (_num_nu > 0)
            {   // synchrotron
                std::ofstream file(_output_sync_spectrum, std::ios_base::trunc);
                file << "# time evolution of synchrotron emission" << '\n';
                file << "# first  line: #num of frequencies, numin, numax, dln" << '\n';
                file << "# frequencies are const in units of critical frequecny which changes with timestep" << '\n';
                file << "# second line: initial time and j(nu)" << '\n';
                file << "# ...2+nth line: nth time and corresponding j(nu)" << '\n';
                file << "# nfrequencies=" << _num_nu << "   numin=" << _numin << "   numax=" << _numax
                     << "   dlnu=" << _dlnu << '\n';

                file.close();
            }
            std::cout << "  output(end): " << '\n';
        }
        else if (a_t > _t_output)
        {
            std::cout << "  output(beg): " << '\n';

            _t_output += _dt_output;
            {
                std::ofstream file(_output_filename, std::ios_base::app);

                file << a_step << "   " << a_t << "   ";
                for (const auto &fp : _fc)
                    file << fp.p << "   ";
                for (const auto &fp : _fc)
                    file << fp.f << "   ";
                file << '\n';

                file.close();
            }
            if (_num_nu > 0)
            {
                std::cout << "    radiation(beg)... ";
                std::ofstream file(_output_sync_spectrum, std::ios_base::app);
                // map uniformly the lagrangian distribution function: prevent roundoff error
                const real_t pmin = _fc.front().p;
                const real_t pmax = _fc.back().p;
                const size_t nodes = _num_nodes / 2;
                const real_t dlp = log(pmax / pmin) / (nodes - 1);
                std::vector<DistFunctNode> fs(nodes, DistFunctNode{zero, zero});
                {
                    auto t_beg = std::chrono::high_resolution_clock::now();
                    fs.front() = _fc.front();
                    fs.back() = _fc.back();
                    size_t j = 0;
                    for (size_t i = 1; i < nodes - 1; ++i)
                    {
                        fs[i].p = pmin * exp(i * dlp);
                        while (_fc[j].p < fs[i].p)
                            ++j;
                        const real_t dfdp = (_fc[j].f - _fc[j - 1].f) / (_fc[j].p - _fc[j - 1].p);
                        fs[i].f = _fc[j - 1].f + dfdp * (fs[i].p - _fc[j - 1].p);
                    }
                    auto t_end = std::chrono::high_resolution_clock::now();
                    _timings["sync_binning"] += t_end - t_beg;
                    //plot_histo(f,100," distribution for synch emission");
                }

                //
                const real_t w = 1.e-12;
                auto fp = [&f = std::as_const(fs), pmin = pmin * (one - w), pmax = pmax * (one + w), dlp](const real_t _p) {
                    if (_p > pmin && _p < pmax)
                    {
                        const size_t i = static_cast<size_t>(floor(log(_p / pmin) / dlp));
                        //assert(i+1<f.size());
                        const real_t dfdp = (f[i + 1].f - f[i].f) / (f[i + 1].p - f[i].p);
                        //assert((dfdp+one)/(dfdp+one)==one);
                        const real_t rv = dfdp * (_p - f[i].p) + f[i].f;
                        //assert(rv/rv==one);
                        return rv;
                    }
                    else if (_p > pmin)
                        return (f.back().f * exp(-_p / pmax));
                    else
                        return (f.front().f * exp(-pmin / _p));
                };

                // get magnetic field is in G
                const real_t B = _icm.get(a_t, "mag-field");
                // critical frequency in GHz
                using namespace PhysConstCgs;
                const real_t nc = 1.e-9 * 1.5 * e * B / (two * Pi * m_e * c);
                // radiation object
                Radiation::Synchrotron rad(fs.front().p, fs.back().p, fp);
                // emission and frequency
                std::vector<real_t> jsy(_num_nu, zero);
                // scale emission freuency around some peak energy
                std::vector<real_t> nsy(_nu.begin(), _nu.end());
                {   // define peak-p at peak of p^4*fs
                    real_t ppk = zero;
                    real_t vpk = zero;
                    for (auto f : fs)
                    {
                        real_t v = pow(f.p, 4) * f.f;
                        if (v > vpk)
                        {
                            vpk = v;
                            ppk = f.p;
                        }
                    }
                    for (auto &n : nsy)
                        n *= pow(ppk, 2);
                }

                {   // thread this time consuming spectrum calculation
                    auto t_beg = std::chrono::high_resolution_clock::now();
                    ThreadUt<const Synchrotron, std::vector<real_t>, const std::vector<real_t>> th_ut{_num_threads,
                                                                                                      _num_nu};
                    th_ut.launch_threads(rad, jsy, nsy);
                    auto t_end = std::chrono::high_resolution_clock::now();
                    _timings["synchrotron"] += t_end - t_beg;
                }
                file << a_step << "   " << a_t << "   ";
                for (auto n : nsy)
                    file << n * nc << "   ";
                for (const auto j : jsy)
                    file << j << "   ";

                file << '\n';

                file.close();
                std::cout << " ...(end): " << '\n';
            }
            std::cout << "  output(end): " << '\n';
        }
    }
} // namespace cr_transport

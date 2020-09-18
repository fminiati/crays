//
// integrator for cr-transport in momentum space
// F. Miniati, (ZH) 08.12.2018
//
#include <vector>
#include <iomanip>
#include <cassert>
#include <map>
#include "Macros.h"

namespace cr_transport
{
    struct CRTranspOp
    {
        CRTranspOp() = default;
        CRTranspOp(CRTranspOp const &) = default;
        CRTranspOp(const real_t n, const real_t T, const real_t B,
                   const real_t dv, const real_t z,
                   const real_t l, const real_t p,
                   const real_t e, const real_t tm,
                   const real_t ut, const real_t cfl)
            : _n(n), _T(T), _B(B), _divv(dv), _zeta(z),
              _lmfp(l), _pmfp(p), _eta(e), _tau_mfp(tm),
              _utime(ut), _cfl(cfl)
        {
        }

        // advance solution by dt
        void operator()(std::vector<DistFunctNode> &a_fc,
                        real_t &a_dt,
                        const size_t a_beg,
                        const size_t a_std,
                        const size_t a_size) const;

    protected:
        // icm plasma model
        real_t _n, _T, _B, _divv, _zeta;

        // transport parameters: mfp scales as: lmfp*(p/pmfp)^eta
        real_t _lmfp, _pmfp, _eta, _tau_mfp;

        // unit of time
        real_t _utime, _cfl;
    };

    struct CRTransport
    {
        // read input and initialization
        void setup(const string a_input);

        // solver
        void solve()
        {
            //
            std::cout << " running from t=" << _t_init << " to " << _t_end << endl;

            // time-integration loop
            real_t t = _t_init;
            unsigned step = 0;
            auto t_beg = std::chrono::high_resolution_clock::now();
            while (t < _t_end && step < _max_step)
            {

                const real_t dt = timestep();
                std::cout << "step=" << step << std::endl;

                try
                {
                    throw(advance(t, dt));
                }
                catch (size_t n)
                {
                    if (n < 2)
                    { // throw it a quit if left with no particles
                        std::cout << "-- Can't integrate single node distribution function --\n program ends here!" << endl;
                        return;
                    }
                }

                ++step;
                t += dt;

                output(t, step);
            }
            auto t_end = std::chrono::high_resolution_clock::now();
            _timings["total"] += t_end - t_beg;

            std::cout << " Timing:" << std::endl;
            const real_t ttot = _timings["total"].count();
            std::cout << "   total:" << ttot << " s" << std::endl;
            for (auto t : _timings)
            {
                if (t.first != "total")
                {
                    const real_t ts = t.second.count();
                    std::cout << "   " << t.first << ":" << ts << "s, or " << ts / ttot << "% of total" << std::endl;
                }
            }
        }

    private:
        // advance solution by dt
        size_t advance(const real_t a_t, const real_t a_dt);

        // manage injection
        void injection(std::vector<DistFunctNode> &a_f, const real_t a_dt);

        // return current timestep
        real_t timestep() const
        {
            assert(_new_dt > zero);
            return _new_dt;
        }

        // write data to output file
        void output(const real_t a_t, const unsigned step);

    private:
        // timing measurements
        std::map<string, std::chrono::duration<real_t>> _timings;

        // number of threads
        size_t _num_threads;

        // distribution function
        std::vector<DistFunctNode> _fc;

        // injection function and time window
        GammaInjection _inj;
        real_t _ti_inj, _te_inj;
        size_t _inj_multiplicity;

        // icm model
        PlasmaModel _icm;

        // timing
        real_t _cfl, _new_dt, _t_init, _t_end, _unit_time_sec;
        size_t _max_step;
        // momentum grid and p boundaries
        size_t _num_nodes;
        real_t _pmin, _pmax;

        // transport parameters: mfp scales as: lmfp*(p/pmfp)^eta
        real_t _lmfp, _pmfp, _eta, _tau_mfp;

        // output file
        string _output_filename, _output_sync_spectrum;
        size_t _num_nu;
        real_t _numin, _numax, _dlnu;
        std::vector<real_t> _nu;
        real_t _dt_output, _t_output;
    };
}; // namespace cr_transport

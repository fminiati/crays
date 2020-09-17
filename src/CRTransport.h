//
// integrator for cr-transport in momentum space
// F. Miniati, (ZH) 08.12.2018
//
#include <vector>
#include <iomanip>
#include <cassert>
#include <forward_list>
#include <map>
#include <random>
#include "mymacros.H"
#include "NewtonRaphson.H"

//
using std::vector;
using std::forward_list;
using std::map;

template <typename C>
void plot_histo(const C& c, const size_t n, const string s=" ") {

  std::cout << s<< endl;
  real_t maxv{};
  for (auto v : c) maxv=std::max(maxv,v.f);
  for (auto v : c) {
    const int x=static_cast<int>(n*v.f/maxv);
    if (x>=1) std::cout << std::fixed
                        << std::setfill('-')
                        << std::setprecision(2)
                        << v.p << " " << string(x,'*') << endl;
  }
};


namespace cr_transport {

  template <typename Op, typename Obj, typename V>
  struct ThreadUt {
    size_t _num_threads;
    size_t _loop_size;

    ThreadUt()=default;
    ThreadUt(ThreadUt const&)=default;
    
    // template thread arg struct
    struct ThreadArg {
      ThreadArg()=default;
      ThreadArg(ThreadArg const&)=default;

      size_t beg, strd, size;
      Op*  op;
      Obj* obj;
      V*   var;
    };

    /// aliases
    using rw_obj_t = std::reference_wrapper<Obj>;
    using rw_var_t = std::reference_wrapper<V>;

    //
    size_t launch_threads(const Op& a_op, Obj& a_obj, vector<V>& a_v) {

      vector<rw_var_t>   v(a_v.begin(),a_v.end());
      vector<rw_obj_t> obj(_num_threads,std::ref(a_obj));
      return launch_threads(a_op,obj,v);
    }

    ///
    size_t launch_threads(const Op& a_op, vector<Obj>& a_obj, V& a_v) {

      vector<rw_var_t> v(_num_threads,std::ref(a_v));
      vector<rw_obj_t> obj(a_obj.begin(),a_obj.end());
      return launch_threads(a_op,obj,v);
    }

    ///
    size_t launch_threads(const Op& a_op, Obj& a_obj, V& a_v) {

      vector<rw_var_t> v(_num_threads,std::ref(a_v));
      vector<rw_obj_t> obj(_num_threads,std::ref(a_obj));
      return launch_threads(a_op,obj,v);
    }

    ///
    size_t launch_threads(const Op& a_op, vector<rw_obj_t>& a_obj, vector<rw_var_t>& a_v) {

      pthread_attr_t th_attr;
      pthread_attr_init(&th_attr);
      pthread_attr_setdetachstate(&th_attr,PTHREAD_CREATE_JOINABLE);

      // create threads
      std::vector<pthread_t> th(_num_threads);
      std::vector<ThreadArg> ta(_num_threads);

      auto t_lambda=[](void* args) ->void* {
                      ThreadArg* ta=static_cast<ThreadArg*>(args);
                      ta->op->operator()(*ta->obj,*ta->var,ta->beg,ta->strd,ta->size);
                      return nullptr;
                    };
      
      size_t err=0;
      for (size_t t=0; t<_num_threads; ++t) {
        ta[t].beg=t;
        ta[t].strd=_num_threads;
        ta[t].size=_loop_size;
        ta[t].op =new Op(a_op);
        ta[t].obj=&a_obj[t].get();
        ta[t].var=&a_v[t].get();
        err += pthread_create(&th[t],&th_attr,t_lambda,(void*)&ta[t]);
      }

      void* status;
      for (auto& t : th) {
        err += pthread_join(t,&status);
      }
      for (auto& a : ta) {
        delete a.op;
        a.op=nullptr;
      }
      assert(err==0);
      return err;
    }
  };

  
  // time dependent model of plasma:
  // n=number density, T=temperature, B=magnetic field
  // divv=average divergence of turbulent velocity on scale _ell
  // zeta=power-law index of dv scaling with separation: dv(l)~dv(l/_ell)^z
  // t=time
  struct PlasmaModel {

    PlasmaModel()
      :
      _it(0) {}

    ~PlasmaModel() {}

    // data
    real_t _ell;
    std::vector<real_t> _n, _T, _B, _divv, _zeta, _t;

    // return variable at time t
    real_t get(const real_t a_t, const string a_var) {
      assert(a_t>=_t.front() && a_t<=_t.back());

      update_index(a_t);

      auto t_interp=[i=_it,t=_t](const auto ta, const auto& f) {
                      return ( (f[i]-f[i-1])/(t[i]-t[i-1]) *(ta-t[i]) + f[i] );
                    };

      if      (a_var=="density")
        return t_interp(a_t,_n);
      else if (a_var=="temperature")
        return t_interp(a_t,_T);
      else if (a_var=="mag-field")
        return t_interp(a_t,_B);
      else if (a_var=="divv-turb")
        return t_interp(a_t,_divv);
      else if (a_var=="zeta-turb")
        return t_interp(a_t,_zeta);

      return zero;
    }

  private:
    // iterpolation index
    void update_index(const real_t a_t) {
      assert(_t.size()>0);

      if ( !(_it>0 && a_t>_t[_it-1] && a_t<_t[_it]) )
        while (a_t >= _t[_it]) ++_it;
    }

    private:
    size_t _it;
  };

  // node for lagrangian description of phase space distribution
  // function
  struct DistFunctNode {
    DistFunctNode()=default;
    DistFunctNode(DistFunctNode const&)=default;
    //DistFunctNode(DistFunctNode &&)=default;

    real_t f;
    real_t p;
  };

  struct Injection {
    Injection()=default;
    Injection(Injection const&)=default;
    virtual ~Injection() {}

    // actual function
    virtual real_t operator () (const real_t x) const = 0;
    // low and high momentum cut-off
    real_t _plo, _phi;
  };

  // f(p)=1/(pi*Gamma) / (1+(p-p0)^2/Gamma^n), normalise to
  // injection time
  struct GammaInjection : public Injection {
    GammaInjection()=default;
    GammaInjection(GammaInjection const&)=default;
    virtual ~GammaInjection() {}

    // actual function
    virtual real_t operator () (const real_t p) const {
      return one/(Pi*_Gamma*_dt) /(one+pow((p-_p0)/_Gamma,_idx));
    }

    size_t _idx;
    real_t _p0, _Gamma, _dt;
  };

  template <typename I>
  struct InjectionOp {
    InjectionOp()=default;
    InjectionOp(const real_t a_pmin,
                const real_t a_dlp,
                const I& a_i)
      :
      _pmin(a_pmin), _dlp(a_dlp), _inject(a_i)
    {}

    //
    void operator() (vector<DistFunctNode>& a_j,
                     const real_t a_dt, const size_t a_beg,
                     const size_t a_stride, const size_t a_size) {

      for (size_t b=a_beg; b<a_size; b+=a_stride) {

        const real_t p=_pmin* exp(b*_dlp);
        a_j[b].p=p;
        if (p>=_inject._plo && p<=_inject._phi)
          a_j[b].f=a_dt*_inject(p);
      }
    }

  protected:
    // p-grid quantities
    real_t _pmin,_dlp;
    // injection funtion: returns f(p)
    I _inject;
  };

  //
  struct CRTranspOp {

    CRTranspOp()=default;
    CRTranspOp(CRTranspOp const&)=default;
    CRTranspOp(const real_t n,  const real_t T, const real_t B,
               const real_t dv, const real_t z, 
               const real_t l, const real_t p,
               const real_t e, const real_t tm,
               const real_t ut, const real_t cfl)
      :
      _n(n), _T(T), _B(B), _divv(dv), _zeta(z),
      _lmfp(l), _pmfp(p), _eta(e), _tau_mfp(tm),
      _utime(ut), _cfl(cfl)
    {}

    // advance solution by dt
    void operator() (vector<DistFunctNode>& a_fc,
                     real_t& a_dt,
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


  struct CRTransport {

    // read input and initialization
    void setup(const string a_input);

    // solver
    void solve() {
      // 
      std::cout << " running from t=" << _t_init << " to " << _t_end <<endl;

      // time-integration loop
      real_t t=_t_init;
      unsigned step=0;
      auto t_beg=std::chrono::high_resolution_clock::now();
      while (t < _t_end && step<_max_step) {

        const real_t dt=timestep();
        std::cout << "step=" << step << std::endl;

        try {
          throw ( advance(t,dt) );
        }
        catch(size_t n) {
          if (n<2) { // throw it a quit if left with no particles
            std::cout << "-- Can't integrate single node distribution function --\n program ends here!" <<endl;
            return;
          }
        }

        ++step;
        t+=dt;

        output(t,step);
      }
      auto t_end=std::chrono::high_resolution_clock::now();
      _timings["total"] += t_end-t_beg;

      std::cout << " Timing:" << std::endl;
      const real_t ttot=_timings["total"].count();
      std::cout << "   total:" << ttot << " s" << std::endl;
      for (auto t : _timings) {
        if (t.first!="total") {
          const real_t ts=t.second.count();
          std::cout << "   "<<t.first<<":"<<ts<<"s, or "<<ts/ttot<<"% of total" <<std::endl;
        }
      }
    }

  private:

    // advance solution by dt
    size_t advance(const real_t a_t, const real_t a_dt);

    // manage injection
    void injection(vector<DistFunctNode>& a_f, const real_t a_dt);

    // return current timestep
    real_t timestep() const {
      assert(_new_dt>zero);
      return _new_dt;
    }

    // write data to output file
    void output (const real_t a_t, const unsigned step);
    
  private:
    // timing measurements
    map<string,std::chrono::duration<real_t>> _timings;

    // number of threads
    size_t _num_threads;

    // distribution function
    vector<DistFunctNode> _fc;

    // injection function and time window
    GammaInjection _inj;
    real_t _ti_inj, _te_inj;
    size_t _inj_multiplicity;

    // icm model
    PlasmaModel _icm;

    // timing
    real_t _cfl,_new_dt,_t_init, _t_end, _unit_time_sec;
    size_t _max_step;
    // momentum grid and p boundaries
    size_t _num_nodes;
    real_t _pmin, _pmax;

    // transport parameters: mfp scales as: lmfp*(p/pmfp)^eta
    real_t _lmfp,_pmfp,_eta,_tau_mfp;

    // output file
    string _output_filename, _output_sync_spectrum;
    size_t _num_nu;
    real_t _numin,_numax,_dlnu;
    vector<real_t> _nu;
    real_t _dt_output, _t_output;
  };
};


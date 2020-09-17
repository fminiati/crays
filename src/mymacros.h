
// somple macros F. Miniati 25.08.06

#ifndef _MY_MACROS_
#define _MY_MACROS_

// my variables
using real_t = double;

template <typename T>
int SGN(const T t) {return (t>T(0) ? 1 : -1);};
template <typename T>
T ABS(const T t) {return (t>T(0) ? t : -t);}


static constexpr auto zero   = (0.0e0);
static constexpr auto half   = (0.5e0);
static constexpr auto one    = (1.0e0);
static constexpr auto two    = (2.0e0);
static constexpr auto three  = (3.0e0);
static constexpr auto four   = (4.0e0);
static constexpr auto five   = (5.0e0);
static constexpr auto six    = (6.0e0);
static constexpr auto seven  = (7.0e0);
static constexpr auto eight  = (8.0e0);
static constexpr auto nine   = (9.0e0);
static constexpr auto ten    = (10.0e0);
static constexpr auto twenty = (20.0e0);
static constexpr auto tenth  = (0.100e0);
static constexpr auto eighth = (0.125e0);
static constexpr auto fifth  = (0.200e0);
static constexpr auto fourth = (0.250e0);
static constexpr auto third  = (one/three);
static constexpr auto Pi     = (3.14159265358979323846e0);
static constexpr auto small  = (1.e-6);
static constexpr auto tiny   = (1.e-9);
static constexpr auto hundred= (1.0e2);
static constexpr auto huge   = (1.0e100);

static constexpr auto MAXITER  = (30);
static constexpr auto TOLERANCE=(1.e-6);
static constexpr auto MAX_EXP  = (20);

#endif

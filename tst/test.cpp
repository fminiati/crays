
#include <cstdio>
#include <functinal>
#include "Romberg.H"
#include "mymacros.H"


int main(int argc, char* argv[]) {

  const int a=4;
  auto f=[&a=a](const real_t x) {
	return a*x*x*x;
	};

   Romberg r;
   std::cout << " int_0^1 x^3 =" << r.integral(zero,one,1.e-9,f) <<std::endl;
}

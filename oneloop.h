//#include "classes.h"
#include "s2lseinline.h"
#include "B0.h"
#include "C0.h"
#include "D0.h"
#include "Cplx.h"
#include "li.h"







inline double B0div(double eps) //maybe these are insufficient for NNLO? Need O(eps,eps^2)?
{
   return 1 / eps;
}

inline double A0div(double ms)
{
   return ms;
}

inline double A0fin(double ms)
{
   return (ms * (1 - log(ms)));
}

inline double B0Refin(Double ps, Double m1s, Double m2s) //function B0A defined in B0.cc is actually B0fin-2
{
   return real(B0A(ps, m1s, m2s)) + 2 ;
}

inline double B0Im(Double ps, Double m1s, Double m2s)
{
   return imag(B0A(ps, m1s, m2s) );
}

inline double C0Re(double p10, double p12, double p20, double m0, double m1, double m2)
{
   return real(C0Q(p10, p12, p20, m0, m1, m2));
}

inline double C0Im(double p10, double p12, double p20, double m0, double m1, double m2)
{
   return imag(C0Q(p10, p12, p20, m0, m1, m2));
}

inline double Power(double base, double exponent)
{
   return pow(base, exponent);
}
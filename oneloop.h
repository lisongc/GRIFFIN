/* oneloop.h: header file for B0.cc, C0.cc, D0.cc */

#include "s2lseinline.h"

#define ASYMP_LIMIT_B0N1 1e5
#define ASYMP_LIMIT_B0S 1e5
#define ASYMP_LIMIT_DeltaB0N1 1e10
#define CPLX_LIMIT_F 1e-8
#define ASYMP_LIMIT_DB01 1e-9
#define SERIES_LIMIT 1e15

// one-loop scalar self-energy,
// B0A = B0(ps,m1s,m2s) -1/delta +EulerGamma -log(4 Pi mu^2) -2
Cplx B0A(Double ps, Double m1s, Double m2s);
Cplx B0A(Cplx ps, Double m1s, Double m2s);

// with one vanishing mass:
Cplx B0MZ(Double ps, Double ms);

// subtracted B0(ps,m1s,m2s)-B0(ps,m1s,0):
Cplx B0N1(Double ps, Double m1s, Double m2s);

// the derivative of B0 with respect to the momentum
Cplx DB0(Double ps, Double m1s, Double m2s);
Cplx DDB0(Double ps, Double m1s, Double m2s);

// the derivative of B0 with respect to the first (squared) mass
Cplx DM1B0(Double ps, Double m1s, Double m2s);

// the derivative of B0 with respect to the first (squared) mass
// and the momentum
Cplx DDM1B0(Double ps, Double m1s, Double m2s);

// some representations of B0 and B0D with subtracted asymptotic for
// large m2s:
//
// B0B = B0 - Delta - log(mu^2) -1 + log(m2^2) 
Cplx B0B(Double ps, Double m1s, Double m2s);

// B0C = B0 - Delta - log(u^2) -1 + log(m2^2) - (m1^2/m2^2) log(m1^2/m2^2)
Cplx B0C(Double ps, Double m1s, Double m2s);

// B0D = B0 - Delta - log(mu^2) -1 + log(m2^2) - (m1^2/m2^2) log(m1^2/m2^2)
//        - ps/(2*m2s)
Cplx B0D(Double ps, Double m1s, Double m2s);

// the derivative of B0 with respect to the momentum,
// with 1/(2 m2^2) subtracted
Cplx DB0D(Double ps, Double m1s, Double m2s);

inline double A0fin(double ms)
{ if(ms==0.) 
    return 0.;
  else
    return(ms*(1-log(ms))); 
}

inline Cplx B0fin(double ps, double m1s, double m2s)
{ return(B0A(ps,m1s,m2s) +2.); }

// the one-loop scalar triangle function
Cplx C0(Double p10, Double p20, Double p12, Double m1, Double m2, Double m3);

// the derivative of C0 with respect to the third (squared) mass
Cplx DM3C0(Double p10, Double p20, Double p12, Double m1, Double m2, Double m3);

//the one-loop scalar box function
//implementation according to A. Denner, Fortschr. Phys. 41(1993) 307
Cplx D0(double ps10, double ps20, double ps30, double ps40, double s12, double s23, double msq0, double msq1, double msq2, double msq3);

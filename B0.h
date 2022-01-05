/* ------------------------------------------------------------- */
/* B0.h, B0.cc                                                   */
/* Stefan Bauberger (stefan@bauberger.de)                        */
/* last revision 4.12.95                                         */
/* ------------------------------------------------------------- */
/* scalar one-loop function B0, some derivatives and the         */
/* discontinuity, for use with the package s2lse                 */
/* ------------------------------------------------------------- */

// discontinuity of B0(s,m1se,ms2s), divided by 2*Pi*I:
// Double DeltaB0OnTwoPiI(Double s, Double m1se, Double m2se)

// derivative of the discontinuity of B0, divided by 2*Pi*I:
// Double DDeltaB0OnTwoPiI(Double s, Double m1se, Double m2se)

// subtracted discontinuity DeltaB0(..)-DeltaB0(m1s=0), divided by 2*Pi*I:
//
// Double DeltaB0OnTwoPiIN1(Double s, Double m1se, Double m2se)

// subtracted discontinuity DeltaB0(..)-DeltaB0(m1=m2=0), divided by 2*Pi*I:
//
// Double DeltaB0OnTwoPiIN2(Double s, Double m1se, Double m2se)

// one-loop scalar self-energy,
// B0A = B0(ps,m1s,m2s) -1/delta +EulerGamma -log(4 Pi mu^2) -2
//
// Cplx B0A(Double ps, Double m1s, Double m2s);

// the same for complex momentum:
//
// Cplx B0A(Cplx ps, Double m1s, Double m2s);

// with one vanishing mass:
//
// Cplx B0MZ(Double ps, Double ms);

// subtracted B0(ps,m1s,m2s)-B0(ps,m1s,0):
//
// Cplx B0N1(Double ps, Double m1s, Double m2s);

// the derivative of B0 with respect to the momentum
//
// Cplx DB0(Double ps, Double m1s, Double m2s);

// the derivative of B0 with respect to the first (squared) mass
//
// Cplx DM1B0(Double ps, Double m1s, Double m2s)

// the derivative of B0 with respect to the first (squared) mass
// and the momentum
// 
// Cplx DDM1B0(Double ps, Double m1s, Double m2s);

// some representations of B0 and B0D with subtracted asymptotic for
// large m2s:
//
// B0B = B0 - Delta - log(mu^2) -1 + log(m2^2) 
// Cplx B0B(Double ps, Double m1s, Double m2s);
// B0C = B0 - Delta - log(u^2) -1 + log(m2^2) - (m1^2/m2^2) log(m1^2/m2^2)
// Cplx B0C(Double ps, Double m1s, Double m2s);
// B0D = B0 - Delta - log(mu^2) -1 + log(m2^2) - (m1^2/m2^2) log(m1^2/m2^2)
//        - ps/(2*m2s)
// Cplx B0D(Double ps, Double m1s, Double m2s);
// the derivative of B0 with respect to the momentum,
// with 1/(2 m2^2) subtracted
// Cplx DB0D(Double ps, Double m1s, Double m2s);

#include "s2lseinline.h"

#ifdef QUADPRECISION
static doubledouble ASYMP_LIMIT_B0N1=SDouble("1e15");
static doubledouble ASYMP_LIMIT_B0S=SDouble("1e15");
static doubledouble ASYMP_LIMIT_DeltaB0N1=SDouble("1e20");
static doubledouble CPLX_LIMIT_F=SDouble("1e-18");
static doubledouble ASYMP_LIMIT_DB01=SDouble("1e-22");
static doubledouble SERIES_LIMIT=SDouble("1e30");
#else
#define ASYMP_LIMIT_B0N1 1e5
#define ASYMP_LIMIT_B0S 1e5
#define ASYMP_LIMIT_DeltaB0N1 1e10
#define CPLX_LIMIT_F 1e-8
#define ASYMP_LIMIT_DB01 1e-9
#define SERIES_LIMIT 1e15
#endif


// discontinuity of B0, divided by 2 Pi I
inline Double DeltaB0OnTwoPiI(Double s, Double m1se, Double m2se)
{ return(sla(s,m1se,m2se)/s); }


// derivative of the discontinuity of B0, divided by 2 Pi I
inline Double DDeltaB0OnTwoPiI(Double s, Double m1se, Double m2se)
{
  return((s*(m1se+m2se)-sqr(m1se-m2se))/(sqr(s)*sla(s,m1se,m2se)));
}

// subtracted discontinuity DeltaB0(..)-DeltaB0(m1s=0), divided by 2 Pi I
inline Double DeltaB0OnTwoPiIN1(Double s, Double m1se, Double m2se)
{
  if (s>ASYMP_LIMIT_DeltaB0N1*(m1se+m2se)) return(-m1se/s);
  else return((m2se-s + slatheta(s,m1se,m2se))/s);
}

// subtracted discontinuity DeltaB0(..)-DeltaB0(m1=m2=0), divided by 2*pi*i
inline Double DeltaB0OnTwoPiIN2(Double s, Double m1se, Double m2se)
{
  if (s>ASYMP_LIMIT_DeltaB0N1*(m1se+m2se)) return(-(m1se+m2se)/s);
  return((slatheta(s,m1se,m2se)-s)/s);
}

Cplx B0A(Double ps, Double m1s, Double m2s);
Cplx B0A(Cplx ps, Double m1s, Double m2s);
Cplx B0MZ(Double ps, Double ms);
Cplx B0N1(Double ps, Double m1s, Double m2s);
Cplx DB0(Double ps, Double m1s, Double m2s);
Cplx DDB0(Double ps, Double m1s, Double m2s);
Cplx DM1B0(Double ps, Double m1s, Double m2s);
Cplx DDM1B0(Double ps, Double m1s, Double m2s);
Cplx B0B(Double ps, Double m1s, Double m2s);
Cplx B0C(Double ps, Double m1s, Double m2s);
Cplx B0D(Double ps, Double m1s, Double m2s);
Cplx DB0D(Double ps, Double m1s, Double m2s);

#ifndef __B0__
#define __B0__


#endif





#ifndef __s2lseinline__
#define __s2lseinline__

#include "Cplx.h"

inline Cplx clog1(Double x) {if (x>0) return(log(x)); else return(log(-x)+PiI);}
   // clog1 = log(x + i*eps) for real x
inline Cplx clog2(Double x) {if (x>0) return(log(x)); else return(log(-x)-PiI);}
   // clog2 = log(x - i*eps) for real x

inline Cplx csqrt1(Double x)    // csqrt1 = sqrt(x+i*eps)
{
  if (x<0) return(I*sqrt(-x)); else return(sqrt(x));
}

inline Cplx csqrt2(Double x)    // csqrt2 = sqrt(x-i*eps)
{
  if (x<0) return(-I*sqrt(-x)); else return(sqrt(x));
}

inline Double max(Double x,Double y)
{
  if (x>y) return(x); else return(y);
}

inline Double max(Double x,Double y,Double z)
{
  if (x>y) { if (x>z) return(x); else return(z);}
  else     { if (y>z) return(y); else return(z);};
}

inline Double min(Double x,Double y)
{
  if (x<y) return(x); else return(y);
}

inline void interchange(Double *x,Double *y)
{
  Double z=*x;
  *x=*y;
  *y=z;
}

//#ifndef __CMATH__
#ifndef __Cplx__
inline Double abs(Double x) {return( fabs(x) );}
#endif

inline Cplx Lm(Double ms, Double mus)   // definition according to ref. [4]
{
  return(EULERGAMMA+clog2(ms/(4.*Pi*mus)));
}

inline Cplx Lm(Double ms) // with mus=1/(4*Pi)
{
  return(EULERGAMMA+clog2(ms));
}

inline Cplx Lp(Double ps, Double mus) { return(Lm(-ps,mus)); }

inline Cplx Lp(Double ps) { return(Lm(-ps)); }

inline Double sqr(const Double& x) { return x*x; }

inline Double trip(const Double& x) { return x*x*x; }

/******************************************************************/
/* the Kaellen function and its square root:                      */
/******************************************************************/

inline Double la(Double x, Double y, Double z) {return(sqr(x-y-z)-4*y*z); }
   // la = lambda - function
inline Double sla(Double x, Double y, Double z) {return(sqrt(la(x,y,z)));}
   // sla = sqrt of lambda function
inline Double snla(Double x, Double y, Double z) {return(sqrt(-la(x,y,z)));}
   // snla = sqrt of negative lambda function

inline Double slatheta(Double x, Double y, Double z)
{
  Double l=la(x,y,z);
  if (l<=Double(0)) return(Double(0));
  if (x<y+z) return(Double(0));
  return(sqrt(l));
}


#endif

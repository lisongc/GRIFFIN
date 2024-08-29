/*-----------------------------------------------------------------------------
li.cc
Ayres Freitas (afreitas@pitt.edu)
last revision: 28 Aug 2024
-------------------------------------------------------------------------------
Calculation of di-logarithm Li_2 using series in terms of logarithms, 
see e.g. arXiv:2201.01678 eq.(9)
-----------------------------------------------------------------------------*/

#include "li.h"

namespace griffin {

static const Double bernoullioverkfac[10] = {
  -1/36., 1/3600., -1/211680., 1/10886400., -1/526901760., 
  691/16999766784000., -1/1120863744000., 3617/181400588328960000., 
 -43867/97072790126247936000., 174611/16860010916664115200000.};

// li2(Double x)=li2(x-I*eps)
Cplx li2(Double x)
{
  Double l, l2, res;
  Cplx lc;
  int i;
  if(x == Double(1)) return(PiSonSix);
  if(x > Double(1)) { lc = log(x)+PiI;  return(-li2(1/x) - PiSonSix - lc*lc/2); }
  if(x < -Double(1)) { l = log(-x);  return(-li2(1/x) - PiSonSix - l*l/2); }
  if(x > 0.51) { return(-li2(1-x) + PiSonSix - log(x)*log(1-x)); }
  l = log(1-x);
  l2 = l*l;
  res = -l - l2/4;
  for(i=0; i<10; i++)
  {
    l *= l2;
    res += bernoullioverkfac[i]*l;
  }
  return(res);
}

// li2c(Double x)=li2(x+I*eps)
Cplx li2c(Double x)
{ return(conj(li2(x))); }

Cplx li2(Cplx x)
{
  Cplx l, l2, res;
  int i;
  if(imag(x) == Double(0)) return(li2(real(x)));
  if(abs(x) > 1.01) { l = log(-x);  return(-li2(1/x) - PiSonSix - l*l/2); }
  if(real(x) > 0.51) { return(-li2(1-x) + PiSonSix - log(x)*log(1-x)); }
  l = log(1-x);
  l2 = l*l;
  res = -l - l2/4;
  for(i=0; i<10; i++)
  {
    l *= l2;
    res += bernoullioverkfac[i]*l;
  }
  return(res);
}

} // namespace

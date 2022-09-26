/*-----------------------------------------------------------------------------
tools.cc
Lisong Chen (lic114@pitt.edu), Ayres Freitas (afreitas@pitt.edu)
last revision: 19 Sep 2022
-------------------------------------------------------------------------------
some useful functions for computing decay widths and cross-sections
-----------------------------------------------------------------------------*/

#include "tools.h"
#include "ff0.h"

// Note: partzwidth will modify the FA and FV objects passed to it,
// by setting the fermion type and input
double partzwidth(FA_SMLO& fa, FV_SMLO& fv, const int type, 
                  const inval& input, int scheme)
{
  double alp = input.get(al)/Pi * (1+ input.get(Delal)), 
         alsp = input.get(als)/Pi;
  double x = sqr(input.get(MZ)/input.get(MT));
  double mf, fac, y, RV, RA;
  switch(type)
  {
    case ELE: mf = input.get(ME);  alsp = 0;  fac = 1;  break;
    case MUO: mf = input.get(MM);  alsp = 0;  fac = 1;  break;
    case TAU: mf = input.get(ML);  alsp = 0;  fac = 1;  break;
    case NUE:
    case NUM:
    case NUT: mf = 0;  alsp = 0;  fac = 1;  break;
    case UQU: mf = input.get(MU);  fac = 3;  break;
    case CQU: mf = input.get(MC);  fac = 3;  break;
    case DQU: mf = input.get(MD);  fac = 3;  break;
    case SQU: mf = input.get(MS);  fac = 3;  break;
    case BQU: mf = input.get(MB);  fac = 3;  break;
    default: return(0);
  }
  // QCD+QED corrections from Phys.Rep.277,189; PRL 101,012002; PRL 108,222003;
  //  PLB 287,209
  y = sqr(mf/input.get(MZ));
  RV = 1+ 0.75*sqr(Qf[type])*alp + alsp
       + (-1.1529539789360381 + x*(0.13037037037037036 
           - 0.02962962962962963*log(x)))*sqr(Qf[type])*sqr(alp)
       - 0.25*sqr(Qf[type])*alp*alsp
       + (1.409230409109778 + x*(0.06518518518518518
           - 0.014814814814814815*log(x)))*sqr(alsp)
       - 12.767064735769054*powint(alsp,3) - 80.0075*powint(alsp,4);
  RA = RV - 2*I3f[type]*(
        (-3.083333 + 0.086420*x + 0.013167*x*x + log(x))*sqr(alsp)
	+ (-15.98773 + 3.72222*log(x) + 1.91667*sqr(log(x)))*powint(alsp,3)
	+ (49.0309 - 17.6637*log(x) + 14.6597*sqr(log(x)) 
	   + 3.6736*powint(log(x),3))*powint(alsp,4));
  RV += 12*y*alsp - 6*sqr(y);
  RA += -6*y - 22*y*alsp + 6*sqr(y);
  fa.setftype(type);
  fa.setinput(input);
  fv.setftype(type);
  fv.setinput(input);
  if(scheme == RUNWIDTHSCHEME)
    fac *= sqrt(1+sqr(input.get(GamZ)/input.get(MZ)));
  return(input.get(MZ)*fac/(12*Pi) 
  	 * (RA*real(fa.result()) + RV*real(fv.result())));
}

// Note: zwidth will modify the FA and FV objects passed to it,
// by setting the fermion type and input
double zwidth(FA_SMLO& fa, FV_SMLO& fv, const inval& input, int scheme)
{
  int gen;
  double res=0;
  for(gen=0; gen<=2; gen++)
    res += partzwidth(fa, fv, ELE+2*gen, input, scheme) 
          + partzwidth(fa, fv, NUE+2*gen, input, scheme)
	  + partzwidth(fa, fv, UQU+2*gen, input, scheme) 
          + partzwidth(fa, fv, DQU+2*gen, input, scheme);
  return(res);
}

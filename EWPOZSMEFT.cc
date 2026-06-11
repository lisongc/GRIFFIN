/*-----------------------------------------------------------------------------
EWPOZSMEFT.cc
Lisong Chen (lic114@pitt.edu), Ayres Freitas (afreitas@pitt.edu), Jackson Wallace (ejw108@pitt.edu)
last revision: 17 Feb 2026
-------------------------------------------------------------------------------
classes for F_A and sw_eff form factors with SMEFT LO corrections
-----------------------------------------------------------------------------*/

#include "EWPOZSMEFT.h"
#include "ff0.h"
#include "ff.h"
#include "oneloop.h"

namespace griffin {

Cplx SW_SMEFTLO::resSMEFTLO(void) const
{
  return(-vz0(ftyp,*ival)/(4*fabs(Qf[ftyp])*az0(ftyp,*ival))*(z0SMEFT(ftyp,VEC,*ival)/vz0(ftyp,*ival)-z0SMEFT(ftyp,AXV,*ival)/az0(ftyp,*ival))); // delta sin(theta_{f,eff}) = -1/(4|Q_f|)v_f(0)^Z/a_f(0)^Z(delta v_f/v_f(0)^Z-delta a_f/a_f(0)^Z).
                                                                                                                                                 // z0SMEFT calculates these deltas, which are el/(2*sw*cw)*(delta g_L^{Zf} +/- delta g_R^{Zf})
}

Cplx FA_SMEFTLO::resSMEFTLO(void) const
{
  return(2*fabs(az0(ftyp,*ival))*fabs(z0SMEFT(ftyp,AXV,*ival))); // In the SMEFT, F_A^f = |a_f^Z|^2 ~ |a_f(0)^Z|^2+2|delta a_f||a_f(0)^Z| 
}

Cplx FV_SMEFTLO::result(void) const
{
  double QVf = 1-4*fabs(Qf[ftyp])*realreg(sw->result());
  return(fa->result()*(QVf*QVf));
}

} // namespace

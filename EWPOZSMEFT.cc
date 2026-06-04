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
  return(-vz0(ftyp,*ival)/(4*fabs(Qf[ftyp])*az0(ftyp,*ival))*(fabs(2*az0(ftyp, *ival))*z0SMEFT(ftyp,VEC,*ival)/vz0(ftyp,*ival)-fabs(2*az0(ftyp, *ival))*z0SMEFT(ftyp,AXV,*ival)/az0(ftyp,*ival))); //remove minus sign to align with convention?
}

Cplx FA_SMEFTLO::resSMEFTLO(void) const
{
  return(2*fabs(az0(ftyp,*ival))*fabs(z0SMEFT(ftyp,AXV,*ival)));
}

Cplx FV_SMEFTLO::result(void) const
{
  double QVf = 1-4*fabs(Qf[ftyp])*realreg(sw->result());
  return(fa->result()*(QVf*QVf));
}

} // namespace

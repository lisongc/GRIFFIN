/*-----------------------------------------------------------------------------
xscda2.cc
Ayres Freitas (afreitas@pitt.edu)
last revision: 5 Apr 2025
-------------------------------------------------------------------------------
provide \Delta\alpha^2 correction to photon s-channel exchange amplitude
-----------------------------------------------------------------------------*/

#include "xscda2.h"
#include "ff0.h"

namespace griffin {

Cplx mat_SMda2::resoffZ2f(void) const
{
  double deltaAlpha = ival->get(Delal);
  double gie0 = g0(it,iff,*ival), gjf0 = g0(ot,off,*ival);
  Cplx sa2 = -deltaAlpha*deltaAlpha*s;
  Cplx mats1 = ((-gie0*gjf0*sa2/s)/s);
  return(mats1);
}

Cplx mat_SMda2::result(void) const
{
  return(mat_SMNNLO::result() + resoffZ2f());
}

} // namespace

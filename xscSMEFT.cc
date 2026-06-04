/*-----------------------------------------------------------------------------
xscSMEFT.cc
E. Jackson Wallace (ejw108@pitt.edu)
last revision: 6 May 2026
-------------------------------------------------------------------------------
matrix element class needed for description of cross-section at LO in the SMEFT away from Z pole
-----------------------------------------------------------------------------*/

#include "xscSMEFT.h"
#include "ff0.h"
#include "ff.h"

namespace griffin {

double mat_SMEFTLO::resoff4f(void) const
{
  double fourFmats =  r4fSMEFT(it, ot, iff, off, *ival);
  return(fourFmats);
}

Cplx mat_SMEFTLO::result(void) const
{
  return(mat_SMNNLO::result()+resoff4f());
}

} // namespace

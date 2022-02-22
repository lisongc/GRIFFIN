/*-----------------------------------------------------------------------------
EWPOZ.cc
Lisong Chen (lic114@pitt.edu), Ayres Freitas (afreitas@pitt.edu)
last revision: 10 Feb 2022
-------------------------------------------------------------------------------
classes for F_A and sw_eff form factors with SM NLO corrections
-----------------------------------------------------------------------------*/

#include "EWPOZ.h"
#include "ff0.h"
#include "ff.h"
#include "oneloop.h"

double SW_SMNLO::res1f(void) const
{
  return(-(az0(ftyp,*ival)*rz1f(ftyp,VEC,*ival) 
          - vz0(ftyp,*ival)*rz1f(ftyp,AXV,*ival))/sqr(az0(ftyp,*ival))
	 /(4*fabs(Qf[ftyp]))); 
}

double SW_SMNLO::res1b(void) const
{
  return(-(az0(ftyp,*ival)*rz1b(ftyp,VEC,*ival) 
          - vz0(ftyp,*ival)*rz1b(ftyp,AXV,*ival))/sqr(az0(ftyp,*ival))
	 /(4*fabs(Qf[ftyp]))); 
}

double FA_SMNLO::res1f(void) const
{
  return(2*az0(ftyp,*ival)*rz1f(ftyp,AXV,*ival)
         - sqr(az0(ftyp,*ival))*rsz1fp(*ival));
}

double FA_SMNLO::res1b(void) const
{
  return(2*az0(ftyp,*ival)*rz1b(ftyp,AXV,*ival)
         - sqr(az0(ftyp,*ival))*rsz1bp(*ival));
}

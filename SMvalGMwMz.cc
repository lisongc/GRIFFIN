/*-----------------------------------------------------------------------------
SMvalGMwMz.cc
Ayres Freitas (afreitas@pitt.edu)
last revision: 21 Jul 2025
-------------------------------------------------------------------------------
input class for the Gmu-MW-MZ scheme, which computes alpha from these inputs;
MW and MZ are in the complex-pole mass scheme
-----------------------------------------------------------------------------*/

#include "deltar.h"
#include "SMvalGMwMz.h"

namespace griffin {

void SMvalGMwMz::compute(void)
{
  SMval::compute();
  
  double ALold = 0, ALcalc = 1/137.;
  double MZs = sqr(data[MZc]), MWs = sqr(data[MWc]), GF = data[Gmu];
  dr_SMNNLO dr(*this);
    
  if(isfinite(data[MZ]*data[MH]*data[MT]*data[MB]*data[MW]*data[als]
              *data[Delal]*data[Gmu])) // only proceed if all parameters needed
	                               // for Delta_r are set
  {
    while(fabs(ALcalc-ALold) > 1e-8*ALcalc)  
     				   // demand 1e-6 technical precision for alpha
    {
      ALold = ALcalc;
      data[al] = ALcalc;
      ALcalc = (0.4501581580785531*GF*MWs*(1 - MWs/MZs))/(1+real(dr.result()));
    }
  }
}

} // namespace

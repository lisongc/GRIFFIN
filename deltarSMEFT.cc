/*-----------------------------------------------------------------------------
deltar.cc
Ayres Freitas (afreitas@pitt.edu)
last revision: 16 Feb 2022
-------------------------------------------------------------------------------
classes for the correction (\Delta r) to the Fermi constant in the SM; 
as well as input class that computes M_W from G_Fermi
-----------------------------------------------------------------------------*/

#include "deltarSMEFT.h"

namespace griffin {

double dr_SMEFTLO::resSMEFT(void) const
{
  double Lambdas = sqr(246),
         CLL0 = ival->get(Cll,1,2,2,1),
         CLL1 = ival->get(Cll,2,1,1,2),
         CPHIL30 = ival->get(Cphil3,1,1),
         CPHIL31 = ival->get(Cphil3,2,2),
         CPHIWB = ival->get(CphiWB),
         CPHID = ival->get(CphiD),
         ALPHA = ival->get(al),
         MZs = sqr(ival->get(MZ)),
         GMU = ival->get(Gmu);
  double cos2sin2 = Pi*ALPHA/(sqrt(2)*GMU*MZs);
  return(1/(Lambdas*GMU)*((CLL0+CLL1)/(2*sqrt(2))-(CPHIL30+CPHIL31)*sqrt(2)/2-sqrt(2)*Pi*ALPHA*CPHIWB/(sqrt(sqrt(2)*GMU*MZs*Pi*ALPHA)*(0.5-sqrt(0.25-cos2sin2)))+CPHID*(1-1/cos2sin2*(0.5+sqrt(0.25-cos2sin2)))/(2*sqrt(2))));
}


void invalGmuSMEFT::compute(void)
{
  double MWsold = 0, MWscalc = 80*80;
  double MZs = sqr(data[MZ]), GF = data[Gmu], alpha = data[al];
  dr_SMEFTLO dr(*this);

  if(isfinite(data[MZ]*data[MH]*data[MT]*data[MB]*data[al]*data[als]*data[Delal]*data[Gmu]
    *dataSMEFTG[Cll][0][1][1][0]*dataSMEFTG[Cll][1][0][0][1]*dataSMEFTG[Cphil3][0][0][0][0]*dataSMEFTG[Cphil3][1][1][0][0]*dataSMEFTG[CphiD][0][0][0][0]*dataSMEFTG[CphiWB][0][0][0][0])) // only proceed if all parameters needed
	                                                                                                                         // for Delta_r are set
  {
    while(fabs(MWscalc-MWsold) > 1e-4)  // demand keV technical precision for m_W
    {
      MWsold = MWscalc;
      data[MW] = sqrt(MWscalc);
      MWscalc = MZs * (0.5 + sqrt(0.25 - 2.221441469079183*alpha/(GF*MZs) 
    			  * (1+real(dr.result()))));
    }
  }
}

} // namespace

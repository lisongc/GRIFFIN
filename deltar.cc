/*-----------------------------------------------------------------------------
deltar.cc
Ayres Freitas (afreitas@pitt.edu)
last revision: 16 Feb 2022
-------------------------------------------------------------------------------
classes for the correction (\Delta r) to the Fermi constant in the SM; 
as well as input class that computes M_W from G_Fermi
-----------------------------------------------------------------------------*/

#include "oneloop.h"
#include "linex.h"
#include "deltar.h"

#include "dr.in"
#include "dr3.in"


#include "draas.grid"

double dr_SMNNLO::res2aas(void) const
{
  double MZs = sqr(ival->get(MZ)),
         MWs = sqr(ival->get(MW)),
         MTs = sqr(ival->get(MT)),
         MBs = sqr(ival->get(MB)),
         ALS = ival->get(als);
// compute O(al als) correction via linear interpolation of a numerical grid
  double z = linex3d(axis1draas, sizeof(axis1draas)/sizeof(double), 
                     axis2draas, sizeof(axis2draas)/sizeof(double), 
                     axis3draas, sizeof(axis3draas)/sizeof(double), 
		     MWs/MZs, MTs/MZs, MBs/MZs, &datadraas[0][0][0]);
  return(drho2aas()/DRHOSCHEME + getalphas(MTs,MZs,ALS)*z);
}

double dr_SMNNLO::res3aas2(void) const
{
  double MZs = sqr(ival->get(MZ)),
         MWs = sqr(ival->get(MW)),
         MTs = sqr(ival->get(MT)),
         AL = ival->get(al) * DRHOSCHEME,
         ALS = ival->get(als);
  // expansion formula from hep-ph/9504413:
  double lz = log(MZs/MTs), xw = (MZs-MWs)/MZs;
  return(-3*MWs/(MZs-MWs) * AL/(1-MWs/MZs) * MTs/(16*Pi*MWs) *
          sqr(getalphas(MTs,MZs,ALS)/Pi) * 
	  (-14.594 + MZs/MTs * (-17.224 + 0.08829*lz + 0.4722*lz*lz +
    	      xw * (22.6367 + 1.2527*lz - 0.8519*lz*lz))
	   + sqr(MZs/MTs) * (-7.7781 - 0.07226*lz + 0.004938*lz*lz +
              xw * (21.497 + 0.05794*lz - 0.006584*lz*lz) -
	      xw*xw * 21.0799))	);
}

#include "dr2fb.grid"

double dr_SMNNLO::res2fb(void) const
{
  double MZs = sqr(ival->get(MZ)),
         MWs = sqr(ival->get(MW)),
         MTs = sqr(ival->get(MT)),
         MHs = sqr(ival->get(MH)),
         deltaAlpha = ival->get(Delal);
// compute O(al_f al_b) correction via linear interpolation of a numerical grid
  double r = linex3d(axis1dr2fb, sizeof(axis1dr2fb)/sizeof(double), 
                     axis2dr2fb, sizeof(axis2dr2fb)/sizeof(double), 
                     axis3dr2fb, sizeof(axis3dr2fb)/sizeof(double), 
                     sqr(log10(MHs)/2), MWs/MZs, MTs/MZs, &datadr2fb[0][0][0]);
  return(r*MTs/MZs + 2*deltaAlpha*res1b());
}

#include "dr2bb.grid"

double dr_SMNNLO::res2bb(void) const
{
  double MZs = sqr(ival->get(MZ)),
         MWs = sqr(ival->get(MW)),
         MHs = sqr(ival->get(MH)),
         deltaAlpha = ival->get(Delal);
// compute O(al_b^2) correction via linear interpolation of a numerical grid
  double r = linex2d(axis1dr2bb, sizeof(axis1dr2bb)/sizeof(double), 
                     axis2dr2bb, sizeof(axis2dr2bb)/sizeof(double), 
                     sqr(log10(MHs)/2), MWs/MZs, &datadr2bb[0][0]);
  return(r);
}

#include "dr3asff.grid"

double dr_SMNNLO::res3ffa2as(void) const
{
  double AL = ival->get(al),
	 ALS = ival->get(als),
	 mz = ival->get(MZ),
	 mw = ival->get(MW),
	 mt = ival->get(MT),
	 deltaAlpha = ival->get(Delal),
	 r1, r2;

  r1 = linex2d(axis1dra2asff, sizeof(axis1dra2asff) / sizeof(double),
               axis2dra2asff, sizeof(axis2dra2asff) / sizeof(double),
               mt / mz, mw / mz,
               &datadra2asff[0][0]);
  r2 = linex2d(axis1dra2asff, sizeof(axis1dra2asff) / sizeof(double),
               axis2dra2asff, sizeof(axis2dra2asff) / sizeof(double),
               mt / mz, mw / mz,
               &datadraasdaff[0][0]);
  return (AL*AL*ALS*r1 + AL*ALS*deltaAlpha*r2);
}   


void invalGmu::compute(void)
{
  double MWsold = 0, MWscalc = 80*80;
  double MZs = sqr(data[MZ]), GF = data[Gmu], alpha = data[al];
  dr_SMNNLO dr(*this);
    
  if(isfinite(data[MZ]*data[MH]*data[MT]*data[MB]*data[al]*data[als]
              *data[Delal]*data[Gmu])) // only proceed if all parameters needed
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

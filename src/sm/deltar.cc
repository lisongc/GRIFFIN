/*
  File: deltar.cc
  Purpose: Evaluate SM Delta-r corrections and solve the Gmu W-mass relation.
  Typed-input integration: Jackson Wallace and Lisong Chen
  Last updated: 2026-11-16
  Original GRIFFIN calculations/author records are retained below or in the
  including source file; this migration changes labels, not the expressions.
*/

/*-----------------------------------------------------------------------------
deltar.cc
Ayres Freitas (afreitas@pitt.edu)
last revision: 16 Feb 2022
-------------------------------------------------------------------------------
classes for the correction (\Delta r) to the Fermi constant in the SM; 
as well as input class that computes M_W from G_Fermi
-----------------------------------------------------------------------------*/

#include "griffin/numerics/oneloop.h"
#include "griffin/numerics/linex.h"
#include "griffin/sm/deltar.h"

namespace griffin {

#include "expressions/dr.in"
#include "expressions/dr3.in"


#include "grids/draas.grid"

double dr_SMNNLO::res2aas(void) const
{
  double MZs = sqr(ival->get(InputPar::ZMassComplex)),
         MWs = sqr(ival->get(InputPar::WMassComplex)),
         MTs = sqr(ival->get(InputPar::topMass)),
         MBs = sqr(ival->get(InputPar::bottomMass)),
         ALS = ival->get(InputPar::alphaS);
// compute O(al als) correction via linear interpolation of a numerical grid
  double z = linex3d(axis1draas, sizeof(axis1draas)/sizeof(double), 
                     axis2draas, sizeof(axis2draas)/sizeof(double), 
                     axis3draas, sizeof(axis3draas)/sizeof(double), 
		     MWs/MZs, MTs/MZs, MBs/MZs, &datadraas[0][0][0]);
  return(drho2aas()/DRHOSCHEME + getalphas(MTs,MZs,ALS)*z);
}

double dr_SMNNLO::res3aas2(void) const
{
  double MZs = sqr(ival->get(InputPar::ZMassComplex)),
         MWs = sqr(ival->get(InputPar::WMassComplex)),
         MTs = sqr(ival->get(InputPar::topMass)),
         AL = ival->get(InputPar::alpha) * DRHOSCHEME,
         ALS = ival->get(InputPar::alphaS);
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

#include "grids/dr2fb.grid"

double dr_SMNNLO::res2fb(void) const
{
  double MZs = sqr(ival->get(InputPar::ZMassComplex)),
         MWs = sqr(ival->get(InputPar::WMassComplex)),
         MTs = sqr(ival->get(InputPar::topMass)),
         MHs = sqr(ival->get(InputPar::HiggsMass)),
         deltaAlpha = ival->get(InputPar::DeltaAlpha);
// compute O(al_f al_b) correction via linear interpolation of a numerical grid
  double r = linex3d(axis1dr2fb, sizeof(axis1dr2fb)/sizeof(double), 
                     axis2dr2fb, sizeof(axis2dr2fb)/sizeof(double), 
                     axis3dr2fb, sizeof(axis3dr2fb)/sizeof(double), 
                     sqr(log10(MHs)/2), MWs/MZs, MTs/MZs, &datadr2fb[0][0][0]);
  return(r*MTs/MZs + 2*deltaAlpha*res1b());
}

#include "grids/dr2bb.grid"

double dr_SMNNLO::res2bb(void) const
{
  double MZs = sqr(ival->get(InputPar::ZMassComplex)),
         MWs = sqr(ival->get(InputPar::WMassComplex)),
         MHs = sqr(ival->get(InputPar::HiggsMass)),
         deltaAlpha = ival->get(InputPar::DeltaAlpha);
// compute O(al_b^2) correction via linear interpolation of a numerical grid
  double r = linex2d(axis1dr2bb, sizeof(axis1dr2bb)/sizeof(double), 
                     axis2dr2bb, sizeof(axis2dr2bb)/sizeof(double), 
                     sqr(log10(MHs)/2), MWs/MZs, &datadr2bb[0][0]);
  return(r);
}

#include "grids/dr3asff.grid"

double dr_SMNNLO::res3ffa2as(void) const
{
  double AL = ival->get(InputPar::alpha),
	 ALS = ival->get(InputPar::alphaS),
	 mz = ival->get(InputPar::ZMassComplex),
	 mw = ival->get(InputPar::WMassComplex),
	 mt = ival->get(InputPar::topMass),
	 deltaAlpha = ival->get(InputPar::DeltaAlpha),
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
  double MZs = sqr(data[index(InputPar::ZMassComplex)]), GF = data[index(InputPar::fermiConstant)], alpha = data[index(InputPar::alpha)];
  dr_SMNNLO dr(*this);
    
  if(isfinite(data[index(InputPar::ZMassComplex)]*data[index(InputPar::HiggsMass)]*data[index(InputPar::topMass)]*data[index(InputPar::bottomMass)]*data[index(InputPar::alpha)]*data[index(InputPar::alphaS)]
              *data[index(InputPar::DeltaAlpha)]*data[index(InputPar::fermiConstant)])) // only proceed if all parameters needed
	                               // for Delta_r are set
  {
    while(fabs(MWscalc-MWsold) > 1e-4)  // demand keV technical precision for m_W
    {
      MWsold = MWscalc;
      data[index(InputPar::WMassComplex)] = sqrt(MWscalc);
      MWscalc = MZs * (0.5 + sqrt(0.25 - 2.221441469079183*alpha/(GF*MZs) 
    			  * (1+real(dr.result()))));
    }
  }
}

} // namespace

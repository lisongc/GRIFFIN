/*
  File: EWPOZ2.cc
  Purpose: Evaluate cumulative NNLO+ SM pole form factors using typed inputs.
  Typed-input integration: Jackson Wallace and Lisong Chen
  Last updated: 2026-11-16
  Original GRIFFIN calculations/author records are retained below or in the
  including source file; this migration changes labels, not the expressions.
*/

/*-----------------------------------------------------------------------------
EWPOZ2.cc
Lisong Chen (lic114@pitt.edu), Ayres Freitas (afreitas@pitt.edu)
last revision: 19 Sep 2022
-------------------------------------------------------------------------------
classes for F_A and sw_eff form factors with SM NNLO+ corrections;
-----------------------------------------------------------------------------*/

#include "griffin/sm/EWPOZ2.h"
#include "griffin/sm/ff0.h"
#include "griffin/sm/ff.h"
#include "griffin/numerics/oneloop.h"
#include "griffin/numerics/linex.h"

namespace griffin {

#include "grids/kpaas.grid"

// corrections are computed using linear interpolation of numerical grids

double SW_SMNNLO::res2aas(void) const
{
  double MZs = sqr(ival->get(InputPar::ZMassComplex)),
         MWs = sqr(ival->get(InputPar::WMassComplex)),
         MTs = sqr(ival->get(InputPar::topMass)),
         MBs = sqr(ival->get(InputPar::bottomMass)),
         ALS = ival->get(InputPar::alphaS);
  double z = linex3d(axis1kpaas, sizeof(axis1kpaas)/sizeof(double),
                     axis2kpaas, sizeof(axis2kpaas)/sizeof(double),
                     axis3kpaas, sizeof(axis3kpaas)/sizeof(double),
		     MWs/MZs, MTs/MZs, MBs/MZs, &datakpaas[0][0][0]);
  return(drho2aas()/DRHOSCHEME + getalphas(MTs,MZs,ALS)*(1-MWs/MZs)*z);
  /* write result as a sum of the rho parameter contribution and the remainder,
     so that the remainder is numerically smaller and the interpolation is
     more accurate */
}

#include "expressions/sw.in"
#include "expressions/sw3.in"

#include "grids/kpl2fb.grid"
#include "grids/kpu2fb.grid"
#include "grids/kpd2fb.grid"
#include "grids/kpb2fb.grid"

double SW_SMNNLO::res2fb(void) const
{
  double MZs = sqr(ival->get(InputPar::ZMassComplex)),
         MWs = sqr(ival->get(InputPar::WMassComplex)),
         MTs = sqr(ival->get(InputPar::topMass)),
         MHs = sqr(ival->get(InputPar::HiggsMass)),
         ALS = ival->get(InputPar::alphaS),
         deltaAlpha = ival->get(InputPar::DeltaAlpha);
  double r=0;
  switch(ftyp) {
   case Fermion::electron:
   case Fermion::muon:
   case Fermion::tau:
    r = linex3d(axis1kpl2fb, sizeof(axis1kpl2fb)/sizeof(double),
                axis2kpl2fb, sizeof(axis2kpl2fb)/sizeof(double),
                axis3kpl2fb, sizeof(axis3kpl2fb)/sizeof(double),
                sqr(log10(MHs)/2), MWs/MZs, MTs/MZs, &datakpl2fb[0][0][0]);
    break;
   case Fermion::nuElectron:
   case Fermion::nuMuon:
   case Fermion::nuTau:
    r = NAN;
    break;
   case Fermion::u:
   case Fermion::c:
    r = linex3d(axis1kpu2fb, sizeof(axis1kpu2fb)/sizeof(double),
                axis2kpu2fb, sizeof(axis2kpu2fb)/sizeof(double),
                axis3kpu2fb, sizeof(axis3kpu2fb)/sizeof(double),
                sqr(log10(MHs)/2), MWs/MZs, MTs/MZs, &datakpu2fb[0][0][0]);
    break;
   case Fermion::d:
   case Fermion::s:
    r = linex3d(axis1kpd2fb, sizeof(axis1kpd2fb)/sizeof(double),
                axis2kpd2fb, sizeof(axis2kpd2fb)/sizeof(double),
                axis3kpd2fb, sizeof(axis3kpd2fb)/sizeof(double),
                sqr(log10(MHs)/2), MWs/MZs, MTs/MZs, &datakpd2fb[0][0][0]);
    break;
   case Fermion::b:
    r = linex3d(axis1kpb2fb, sizeof(axis1kpb2fb)/sizeof(double),
                axis2kpb2fb, sizeof(axis2kpb2fb)/sizeof(double),
                axis3kpb2fb, sizeof(axis3kpb2fb)/sizeof(double),
                sqr(log10(MHs)/2), MWs/MZs, MTs/MZs, &datakpb2fb[0][0][0]);
    break;
  }
  return(r*MTs/MZs*(1-MWs/MZs) + deltaAlpha*res1b());
  /* term proportional to \delta\alpha can be treated separately (analytically);
     for the remainder, a prefactor proportional to mt^2 is factored out,
     to reduce the depenendence on mt and thus improve the accuracy of the
     interpolation */
}

#include "grids/kp2bb.grid"

double SW_SMNNLO::res2bb(void) const
{
  double mz = ival->get(InputPar::ZMassComplex),
         mw = ival->get(InputPar::WMassComplex),
         mh = ival->get(InputPar::HiggsMass),
         mt = ival->get(InputPar::topMass),
         deltaAlpha = ival->get(InputPar::DeltaAlpha);
  double r=0;
  switch(ftyp) {
   case Fermion::electron:
   case Fermion::muon:
   case Fermion::tau:  r = linex2d(axis1kpl2bb, sizeof(axis1kpl2bb)/sizeof(double),
	  axis2kpl2bb, sizeof(axis2kpl2bb)/sizeof(double),
	  (mh/mz)/(125.1/91.1876),
			  sqr((mw/mz)/(80.385/91.1876)) -1,
			  &datakpl2bb[0][0]);
              break;
   case Fermion::nuElectron:
   case Fermion::nuMuon:
   case Fermion::nuTau:  r = NAN;
              break;
   case Fermion::u:
   case Fermion::c:  r = linex2d(axis1kpu2bb, sizeof(axis1kpu2bb)/sizeof(double),
	  axis2kpu2bb, sizeof(axis2kpu2bb)/sizeof(double),
	  (mh/mz)/(125.1/91.1876),
			  sqr((mw/mz)/(80.385/91.1876)) -1,
			  &datakpu2bb[0][0]);
              break;
   case Fermion::d:
   case Fermion::s:  r = linex2d(axis1kpd2bb, sizeof(axis1kpd2bb)/sizeof(double),
	  axis2kpd2bb, sizeof(axis2kpd2bb)/sizeof(double),
	  (mh/mz)/(125.1/91.1876),
			  sqr((mw/mz)/(80.385/91.1876)) -1,
			  &datakpd2bb[0][0]);
              break;
   case Fermion::b:  r = linex3d(axis1kpb2bb, sizeof(axis1kpb2bb)/sizeof(double),
	  axis2kpb2bb, sizeof(axis2kpb2bb)/sizeof(double),
	  axis3kpb2bb, sizeof(axis3kpb2bb)/sizeof(double),
	  (mh/mz)/(125.1/91.1876),
			  sqr((mw/mz)/(80.385/91.1876)) -1,
			  sqr((mt/mz)/(173.2/91.1876)) -1,
			  &datakpb2bb[0][0][0]);
              break;
  }
  return(r*(1-sqr(mw/mz)));
}

#include "grids/sw3asff.grid"

double SW_SMNNLO::res3ffa2as(void) const
{
  double AL = ival->get(InputPar::alpha),
	 ALS = ival->get(InputPar::alphaS),
	 mz = ival->get(InputPar::ZMassComplex),
	 mw = ival->get(InputPar::WMassComplex),
	 mt = ival->get(InputPar::topMass),
	 deltaAlpha = ival->get(InputPar::DeltaAlpha),
	 r1, r2;

  r1 = linex2d(axis1swa2asff, sizeof(axis1swa2asff) / sizeof(double),
               axis2swa2asff, sizeof(axis2swa2asff) / sizeof(double),
               mt / mz, mw / mz,
               &dataswa2asff[0][0]);
  r2 = linex2d(axis1swa2asff, sizeof(axis1swa2asff) / sizeof(double),
               axis2swa2asff, sizeof(axis2swa2asff) / sizeof(double),
               mt / mz, mw / mz,
               &dataswaasdaff[0][0]);
  return (AL*AL*ALS*r1 + AL*ALS*deltaAlpha*r2);
}

Cplx SW_SMNNLO::errest(void) const
{
  switch(ftyp) {
   case Fermion::electron:
   case Fermion::muon:
   case Fermion::tau:
   case Fermion::u:
   case Fermion::c:
   case Fermion::d:
   case Fermion::s:  return(4.3e-5);
   case Fermion::b:  return(5.3e-5);   // from 1906.08815
   case Fermion::nuElectron:
   case Fermion::nuMuon:
   case Fermion::nuTau:  return(0);
  }
  return(0);
}


double FA_SMNNLO::delrhofac(void) const
{
  double cws = sqr(ival->get(InputPar::WMassComplex)/ival->get(InputPar::ZMassComplex)),
         alp = ival->get(InputPar::alpha);
  return(alp*Pi*(1-2*cws)/(4*cws*sqr(1-cws)));
}

#include "grids/fa_aas.grid"

double FA_SMNNLO::res2aas(void) const
{
  double MZs = sqr(ival->get(InputPar::ZMassComplex)),
         MWs = sqr(ival->get(InputPar::WMassComplex)),
         MTs = sqr(ival->get(InputPar::topMass)),
         MBs = sqr(ival->get(InputPar::bottomMass)),
         ALS = ival->get(InputPar::alphaS);
  double z = linex3d(axis1fa2aas, sizeof(axis1fa2aas)/sizeof(double),
                     axis2fa2aas, sizeof(axis2fa2aas)/sizeof(double),
                     axis3fa2aas, sizeof(axis3fa2aas)/sizeof(double),
		     MWs/MZs, MTs/MZs, MBs/MZs, &datafa2aas[0][0][0]);
  return(drho2aas()/DRHOSCHEME + getalphas(MTs,MZs,ALS)*z);
  /* write result as a sum of the rho parameter contribution and the remainder,
     so that the remainder is numerically smaller and the interpolation is
     more accurate */
}

#include "expressions/fa.in"
#include "expressions/fa3.in"

#include "grids/fal2fb.grid"
#include "grids/fan2fb.grid"
#include "grids/fau2fb.grid"
#include "grids/fad2fb.grid"
#include "grids/fab2fb.grid"

double FA_SMNNLO::res2fb(void) const
{
  double MZs = sqr(ival->get(InputPar::ZMassComplex)),
         MWs = sqr(ival->get(InputPar::WMassComplex)),
         MTs = sqr(ival->get(InputPar::topMass)),
         MHs = sqr(ival->get(InputPar::HiggsMass)),
         ALS = ival->get(InputPar::alphaS),
         deltaAlpha = ival->get(InputPar::DeltaAlpha);
  double r=0;
  switch(ftyp) {
   case Fermion::electron:
   case Fermion::muon:
   case Fermion::tau:
    r = linex3d(axis1fal2fb, sizeof(axis1fal2fb)/sizeof(double),
                axis2fal2fb, sizeof(axis2fal2fb)/sizeof(double),
                axis3fal2fb, sizeof(axis3fal2fb)/sizeof(double),
                sqr(log10(MHs)/2), MWs/MZs, MTs/MZs, &datafal2fb[0][0][0]);
    break;
   case Fermion::nuElectron:
   case Fermion::nuMuon:
   case Fermion::nuTau:
    r = linex3d(axis1fan2fb, sizeof(axis1fan2fb)/sizeof(double),
                axis2fan2fb, sizeof(axis2fan2fb)/sizeof(double),
                axis3fan2fb, sizeof(axis3fan2fb)/sizeof(double),
                sqr(log10(MHs)/2), MWs/MZs, MTs/MZs, &datafan2fb[0][0][0]);
    break;
   case Fermion::u:
   case Fermion::c:
    r = linex3d(axis1fau2fb, sizeof(axis1fau2fb)/sizeof(double),
                axis2fau2fb, sizeof(axis2fau2fb)/sizeof(double),
                axis3fau2fb, sizeof(axis3fau2fb)/sizeof(double),
                sqr(log10(MHs)/2), MWs/MZs, MTs/MZs, &datafau2fb[0][0][0]);
    break;
   case Fermion::d:
   case Fermion::s:
    r = linex3d(axis1fad2fb, sizeof(axis1fad2fb)/sizeof(double),
                axis2fad2fb, sizeof(axis2fad2fb)/sizeof(double),
                axis3fad2fb, sizeof(axis3fad2fb)/sizeof(double),
                sqr(log10(MHs)/2), MWs/MZs, MTs/MZs, &datafad2fb[0][0][0]);
    break;
   case Fermion::b:
    r = linex3d(axis1fab2fb, sizeof(axis1fab2fb)/sizeof(double),
                axis2fab2fb, sizeof(axis2fab2fb)/sizeof(double),
                axis3fab2fb, sizeof(axis3fab2fb)/sizeof(double),
                sqr(log10(MHs)/2), MWs/MZs, MTs/MZs, &datafab2fb[0][0][0]);
    break;
  }
  return(r*MTs/MZs + 2*deltaAlpha*res1b());
  /* term proportional to \delta\alpha can be treated separately (analytically);
     for the remainder, a prefactor proportional to mt^2 is factored out,
     to reduce the depenendence on mt and thus improve the accuracy of the
     interpolation */
}

#include "grids/fa2bb.grid"

double FA_SMNNLO::res2bb(void) const
{
  double mz = ival->get(InputPar::ZMassComplex),
         mw = ival->get(InputPar::WMassComplex),
         mh = ival->get(InputPar::HiggsMass),
         mt = ival->get(InputPar::topMass),
         deltaAlpha = ival->get(InputPar::DeltaAlpha);
  double r=0;
  switch(ftyp) {
   case Fermion::electron:
   case Fermion::muon:
   case Fermion::tau:  r = linex2d(axis1fal2bb, sizeof(axis1fal2bb)/sizeof(double),
	  axis2fal2bb, sizeof(axis2fal2bb)/sizeof(double),
	  (mh/mz)/(125.1/91.1876),
			  sqr((mw/mz)/(80.385/91.1876)) -1,
			  &datafal2bb[0][0]);
              break;
   case Fermion::nuElectron:
   case Fermion::nuMuon:
   case Fermion::nuTau:  r = linex2d(axis1fan2bb, sizeof(axis1fan2bb)/sizeof(double),
	  axis2fan2bb, sizeof(axis2fan2bb)/sizeof(double),
	  (mh/mz)/(125.1/91.1876),
			  sqr((mw/mz)/(80.385/91.1876)) -1,
			  &datafan2bb[0][0]);
              break;
   case Fermion::u:
   case Fermion::c:  r = linex2d(axis1fau2bb, sizeof(axis1fau2bb)/sizeof(double),
	  axis2fau2bb, sizeof(axis2fau2bb)/sizeof(double),
	  (mh/mz)/(125.1/91.1876),
			  sqr((mw/mz)/(80.385/91.1876)) -1,
			  &datafau2bb[0][0]);
              break;
   case Fermion::d:
   case Fermion::s:  r = linex2d(axis1fad2bb, sizeof(axis1fad2bb)/sizeof(double),
	  axis2fad2bb, sizeof(axis2fad2bb)/sizeof(double),
	  (mh/mz)/(125.1/91.1876),
			  sqr((mw/mz)/(80.385/91.1876)) -1,
			  &datafad2bb[0][0]);
              break;
   case Fermion::b:  r = linex3d(axis1fab2bb, sizeof(axis1fab2bb)/sizeof(double),
	  axis2fab2bb, sizeof(axis2fab2bb)/sizeof(double),
	  axis3fab2bb, sizeof(axis3fab2bb)/sizeof(double),
	  (mh/mz)/(125.1/91.1876),
			  sqr((mw/mz)/(80.385/91.1876)) -1,
			  sqr((mt/mz)/(173.2/91.1876)) -1,
			  &datafab2bb[0][0][0]);
              break;
  }
  return(r);
}

#include "expressions/nfaas.in"

double SW_SMNNLO::res2aasnf(void) const
{
  return(-(az0(ftyp,*ival)*zaas(ftyp,Current::vector,*ival)
          - vz0(ftyp,*ival)*zaas(ftyp,Current::axial,*ival))/sqr(az0(ftyp,*ival))
	 /(4*fabs(Qf[index(ftyp)])));
}

double FA_SMNNLO::res2aasnf(void) const
{
  return(2*az0(ftyp,*ival)*zaas(ftyp,Current::axial,*ival));
}

#include "grids/fa3asff.grid"

double FA_SMNNLO::res3ffa2as(void) const
{
  double AL = ival->get(InputPar::alpha),
	 ALS = ival->get(InputPar::alphaS),
	 mz = ival->get(InputPar::ZMassComplex),
	 mw = ival->get(InputPar::WMassComplex),
	 mt = ival->get(InputPar::topMass),
	 deltaAlpha = ival->get(InputPar::DeltaAlpha),
	 r1, r2;

  r1 = linex2d(axis1fa3a2asff, sizeof(axis1fa3a2asff) / sizeof(double),
               axis2fa3a2asff, sizeof(axis2fa3a2asff) / sizeof(double),
               mt / mz, mw / mz,
               &datafa3a2asff[0][0]);
  r2 = linex2d(axis1fa3a2asff, sizeof(axis1fa3a2asff) / sizeof(double),
               axis2fa3a2asff, sizeof(axis2fa3a2asff) / sizeof(double),
               mt / mz, mw / mz,
               &datafa3adaasff[0][0]);
  return (AL*AL*AL*ALS*r1 +  AL*AL*ALS*deltaAlpha*r2);

}

Cplx FA_SMNNLO::errest(void) const
{
  switch(ftyp) {
   case Fermion::electron:
   case Fermion::muon:
   case Fermion::tau:  return(0.4e-5);
   case Fermion::nuElectron:
   case Fermion::nuMuon:
   case Fermion::nuTau:  return(0.6e-5);
   case Fermion::u:
   case Fermion::c:  return(0.5e-5);
   case Fermion::d:
   case Fermion::s:  return(0.7e-5);
   case Fermion::b:  return(0.3e-5);   // from 1906.08815
  }
  return(0);
}


Cplx FV_SMNNLO::result(void) const
{
  double IVf = (az0(ftyp,*ival)*(iz1f(ftyp,Current::vector,*ival)+iz1b(ftyp,Current::vector,*ival))
	         - vz0(ftyp,*ival)*(iz1f(ftyp,Current::axial,*ival)+iz1b(ftyp,Current::axial,*ival)))
		    /sqr(az0(ftyp,*ival));
  double QVf = 1-4*fabs(Qf[index(ftyp)])*realreg(sw->result());
  return(fa->result()*(QVf*QVf + IVf*IVf));
}

Cplx FV_SMNNLO::errest(void) const
{
  switch(ftyp) {
   case Fermion::electron:
   case Fermion::muon:
   case Fermion::tau:  return(0.5e-5);
   case Fermion::nuElectron:
   case Fermion::nuMuon:
   case Fermion::nuTau:  return(0.6e-5);
   case Fermion::u:
   case Fermion::c:  return(0.5e-5);
   case Fermion::d:
   case Fermion::s:  return(1.1e-5);
   case Fermion::b:  return(1.1e-5);   // from 1906.08815
  }
  return(0);
}

} // namespace

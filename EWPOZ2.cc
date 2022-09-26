/*-----------------------------------------------------------------------------
EWPOZ2.cc
Lisong Chen (lic114@pitt.edu), Ayres Freitas (afreitas@pitt.edu)
last revision: 19 Sep 2022
-------------------------------------------------------------------------------
classes for F_A and sw_eff form factors with SM NNLO+ corrections;
-----------------------------------------------------------------------------*/

#include "EWPOZ2.h"
#include "ff0.h"
#include "ff.h"
#include "oneloop.h"
#include "linex.h"

#include "kpaas.grid"

// corrections are computed using linear interpolation of numerical grids

double SW_SMNNLO::res2aas(void) const
{
  double MZs = sqr(ival->get(MZ)),
         MWs = sqr(ival->get(MW)),
         MTs = sqr(ival->get(MT)),
         MBs = sqr(ival->get(MB)),
         ALS = ival->get(als);
  double z = linex3d(axis1kpaas, sizeof(axis1kpaas)/sizeof(double), 
                     axis2kpaas, sizeof(axis2kpaas)/sizeof(double), 
                     axis3kpaas, sizeof(axis3kpaas)/sizeof(double), 
		     MWs/MZs, MTs/MZs, MBs/MZs, &datakpaas[0][0][0]);
  return(drho2aas()/DRHOSCHEME + getalphas(MTs,MZs,ALS)*(1-MWs/MZs)*z);
  /* write result as a sum of the rho parameter contribution and the remainder,
     so that the remainder is numerically smaller and the interpolation is
     more accurate */
}

#include "sw.in"

#include "kpl2fb.grid"
#include "kpu2fb.grid"
#include "kpd2fb.grid"
#include "kpb2fb.grid"

double SW_SMNNLO::res2fb(void) const
{
  double MZs = sqr(ival->get(MZ)),
         MWs = sqr(ival->get(MW)),
         MTs = sqr(ival->get(MT)),
         MHs = sqr(ival->get(MH)),
         ALS = ival->get(als),
         deltaAlpha = ival->get(Delal);
  double r=0;
  switch(ftyp) {
   case ELE:
   case MUO:
   case TAU:
    r = linex3d(axis1kpl2fb, sizeof(axis1kpl2fb)/sizeof(double), 
                axis2kpl2fb, sizeof(axis2kpl2fb)/sizeof(double), 
                axis3kpl2fb, sizeof(axis3kpl2fb)/sizeof(double), 
                sqr(log10(MHs)/2), MWs/MZs, MTs/MZs, &datakpl2fb[0][0][0]);
    break;
   case NUE:
   case NUM:
   case NUT:
    r = NAN;
    break;
   case UQU:
   case CQU:
    r = linex3d(axis1kpu2fb, sizeof(axis1kpu2fb)/sizeof(double), 
                axis2kpu2fb, sizeof(axis2kpu2fb)/sizeof(double), 
                axis3kpu2fb, sizeof(axis3kpu2fb)/sizeof(double), 
                sqr(log10(MHs)/2), MWs/MZs, MTs/MZs, &datakpu2fb[0][0][0]);
    break;
   case DQU:
   case SQU:
    r = linex3d(axis1kpd2fb, sizeof(axis1kpd2fb)/sizeof(double), 
                axis2kpd2fb, sizeof(axis2kpd2fb)/sizeof(double), 
                axis3kpd2fb, sizeof(axis3kpd2fb)/sizeof(double), 
                sqr(log10(MHs)/2), MWs/MZs, MTs/MZs, &datakpd2fb[0][0][0]);
    break;
   case BQU:
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

#include "kp2bb.grid"

double SW_SMNNLO::res2bb(void) const
{
  double mz = ival->get(MZ),
         mw = ival->get(MW),
         mh = ival->get(MH),
         mt = ival->get(MT),
         deltaAlpha = ival->get(Delal);
  double r=0;
  switch(ftyp) {
   case ELE:
   case MUO:
   case TAU:  r = linex2d(axis1kpl2bb, sizeof(axis1kpl2bb)/sizeof(double), 
                	  axis2kpl2bb, sizeof(axis2kpl2bb)/sizeof(double), 
                	  (mh/mz)/(125.1/91.1876),
			  sqr((mw/mz)/(80.385/91.1876)) -1,
			  &datakpl2bb[0][0]);
              break;
   case NUE:
   case NUM:
   case NUT:  r = NAN;
              break;
   case UQU:
   case CQU:  r = linex2d(axis1kpu2bb, sizeof(axis1kpu2bb)/sizeof(double), 
                	  axis2kpu2bb, sizeof(axis2kpu2bb)/sizeof(double), 
                	  (mh/mz)/(125.1/91.1876),
			  sqr((mw/mz)/(80.385/91.1876)) -1,
			  &datakpu2bb[0][0]);
              break;
   case DQU:
   case SQU:  r = linex2d(axis1kpd2bb, sizeof(axis1kpd2bb)/sizeof(double), 
                	  axis2kpd2bb, sizeof(axis2kpd2bb)/sizeof(double), 
                	  (mh/mz)/(125.1/91.1876),
			  sqr((mw/mz)/(80.385/91.1876)) -1,
			  &datakpd2bb[0][0]);
              break;
   case BQU:  r = linex3d(axis1kpb2bb, sizeof(axis1kpb2bb)/sizeof(double), 
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

Cplx SW_SMNNLO::errest(void) const
{
  switch(ftyp) {
   case ELE:
   case MUO:
   case TAU:
   case UQU:
   case CQU:
   case DQU:
   case SQU:  return(4.3e-5);
   case BQU:  return(5.3e-5);   // from 1906.08815
   case NUE:
   case NUM:
   case NUT:  return(0);
  }
  return(0);
}


double FA_SMNNLO::delrhofac(void) const
{
  double cws = sqr(ival->get(MW)/ival->get(MZ)),
         alp = ival->get(al);
  return(alp*Pi*(1-2*cws)/(4*cws*sqr(1-cws)));
}

#include "fa_aas.grid"

double FA_SMNNLO::res2aas(void) const
{
  double MZs = sqr(ival->get(MZ)),
         MWs = sqr(ival->get(MW)),
         MTs = sqr(ival->get(MT)),
         MBs = sqr(ival->get(MB)),
         ALS = ival->get(als);
  double z = linex3d(axis1fa2aas, sizeof(axis1fa2aas)/sizeof(double), 
                     axis2fa2aas, sizeof(axis2fa2aas)/sizeof(double), 
                     axis3fa2aas, sizeof(axis3fa2aas)/sizeof(double), 
		     MWs/MZs, MTs/MZs, MBs/MZs, &datafa2aas[0][0][0]);
  return(drho2aas()/DRHOSCHEME + getalphas(MTs,MZs,ALS)*z);
  /* write result as a sum of the rho parameter contribution and the remainder,
     so that the remainder is numerically smaller and the interpolation is
     more accurate */
}

#include "fa.in"

#include "fal2fb.grid"
#include "fan2fb.grid"
#include "fau2fb.grid"
#include "fad2fb.grid"
#include "fab2fb.grid"

double FA_SMNNLO::res2fb(void) const
{
  double MZs = sqr(ival->get(MZ)),
         MWs = sqr(ival->get(MW)),
         MTs = sqr(ival->get(MT)),
         MHs = sqr(ival->get(MH)),
         ALS = ival->get(als),
         deltaAlpha = ival->get(Delal);
  double r=0;
  switch(ftyp) {
   case ELE:
   case MUO:
   case TAU:
    r = linex3d(axis1fal2fb, sizeof(axis1fal2fb)/sizeof(double), 
                axis2fal2fb, sizeof(axis2fal2fb)/sizeof(double), 
                axis3fal2fb, sizeof(axis3fal2fb)/sizeof(double), 
                sqr(log10(MHs)/2), MWs/MZs, MTs/MZs, &datafal2fb[0][0][0]);
    break;
   case NUE:
   case NUM:
   case NUT:
    r = linex3d(axis1fan2fb, sizeof(axis1fan2fb)/sizeof(double), 
                axis2fan2fb, sizeof(axis2fan2fb)/sizeof(double), 
                axis3fan2fb, sizeof(axis3fan2fb)/sizeof(double), 
                sqr(log10(MHs)/2), MWs/MZs, MTs/MZs, &datafan2fb[0][0][0]);
    break;
   case UQU:
   case CQU:
    r = linex3d(axis1fau2fb, sizeof(axis1fau2fb)/sizeof(double), 
                axis2fau2fb, sizeof(axis2fau2fb)/sizeof(double), 
                axis3fau2fb, sizeof(axis3fau2fb)/sizeof(double), 
                sqr(log10(MHs)/2), MWs/MZs, MTs/MZs, &datafau2fb[0][0][0]);
    break;
   case DQU:
   case SQU:
    r = linex3d(axis1fad2fb, sizeof(axis1fad2fb)/sizeof(double), 
                axis2fad2fb, sizeof(axis2fad2fb)/sizeof(double), 
                axis3fad2fb, sizeof(axis3fad2fb)/sizeof(double), 
                sqr(log10(MHs)/2), MWs/MZs, MTs/MZs, &datafad2fb[0][0][0]);
    break;
   case BQU:
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

#include "fa2bb.grid"

double FA_SMNNLO::res2bb(void) const
{
  double mz = ival->get(MZ),
         mw = ival->get(MW),
         mh = ival->get(MH),
         mt = ival->get(MT),
         deltaAlpha = ival->get(Delal);
  double r=0;
  switch(ftyp) {
   case ELE:
   case MUO:
   case TAU:  r = linex2d(axis1fal2bb, sizeof(axis1fal2bb)/sizeof(double), 
                	  axis2fal2bb, sizeof(axis2fal2bb)/sizeof(double), 
                	  (mh/mz)/(125.1/91.1876),
			  sqr((mw/mz)/(80.385/91.1876)) -1,
			  &datafal2bb[0][0]);
              break;
   case NUE:
   case NUM:
   case NUT:  r = linex2d(axis1fan2bb, sizeof(axis1fan2bb)/sizeof(double), 
                	  axis2fan2bb, sizeof(axis2fan2bb)/sizeof(double), 
                	  (mh/mz)/(125.1/91.1876),
			  sqr((mw/mz)/(80.385/91.1876)) -1,
			  &datafan2bb[0][0]);
              break;
   case UQU:
   case CQU:  r = linex2d(axis1fau2bb, sizeof(axis1fau2bb)/sizeof(double), 
                	  axis2fau2bb, sizeof(axis2fau2bb)/sizeof(double), 
                	  (mh/mz)/(125.1/91.1876),
			  sqr((mw/mz)/(80.385/91.1876)) -1,
			  &datafau2bb[0][0]);
              break;
   case DQU:
   case SQU:  r = linex2d(axis1fad2bb, sizeof(axis1fad2bb)/sizeof(double), 
                	  axis2fad2bb, sizeof(axis2fad2bb)/sizeof(double), 
                	  (mh/mz)/(125.1/91.1876),
			  sqr((mw/mz)/(80.385/91.1876)) -1,
			  &datafad2bb[0][0]);
              break;
   case BQU:  r = linex3d(axis1fab2bb, sizeof(axis1fab2bb)/sizeof(double), 
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

#include "nfaas.in"

double SW_SMNNLO::res2aasnf(void) const
{
  return(-(az0(ftyp,*ival)*zaas(ftyp,VEC,*ival) 
          - vz0(ftyp,*ival)*zaas(ftyp,AXV,*ival))/sqr(az0(ftyp,*ival))
	 /(4*fabs(Qf[ftyp]))); 
}

double FA_SMNNLO::res2aasnf(void) const
{
  return(2*az0(ftyp,*ival)*zaas(ftyp,AXV,*ival));
}

Cplx FA_SMNNLO::errest(void) const
{
  switch(ftyp) {
   case ELE:
   case MUO:
   case TAU:  return(0.4e-5);
   case NUE:
   case NUM:
   case NUT:  return(0.6e-5);
   case UQU:
   case CQU:  return(0.5e-5);
   case DQU:
   case SQU:  return(0.7e-5);
   case BQU:  return(0.3e-5);   // from 1906.08815
  }
  return(0);
}


Cplx FV_SMNNLO::result(void) const
{
  double IVf = (az0(ftyp,*ival)*(iz1f(ftyp,VEC,*ival)+iz1b(ftyp,VEC,*ival)) 
	         - vz0(ftyp,*ival)*(iz1f(ftyp,AXV,*ival)+iz1b(ftyp,AXV,*ival)))
		    /sqr(az0(ftyp,*ival));
  double QVf = 1-4*fabs(Qf[ftyp])*realreg(sw->result());
  return(fa->result()*(QVf*QVf + IVf*IVf));
}

Cplx FV_SMNNLO::errest(void) const
{
  switch(ftyp) {
   case ELE:
   case MUO:
   case TAU:  return(0.5e-5);
   case NUE:
   case NUM:
   case NUT:  return(0.6e-5);
   case UQU:
   case CQU:  return(0.5e-5);
   case DQU:
   case SQU:  return(1.1e-5);
   case BQU:  return(1.1e-5);   // from 1906.08815
  }
  return(0);
}

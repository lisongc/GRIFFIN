/*-----------------------------------------------------------------------------
xscaas.cc
Ayres Freitas (afreitas@pitt.edu)
last revision: 26 Aug 2024
-------------------------------------------------------------------------------
matrix element class needed for description of cross-section at O(Nf al als) 
precision away from Z pole
-----------------------------------------------------------------------------*/

#include "xscaas.h"
#include "ff0.h"
#include "ff.h"
#include "linex.h"
#include "delrho.h"

namespace griffin {

#include "aasoffz.grid"

// O(aas) corrections to effective Z vertex vector factor 
// (which stem from photon-Z mixing as the only s-dependent contribution)
Cplx zhataas(double s, const inval* ival)
{
  double MZs = sqr(ival->get(MZ)),
         MWs = sqr(ival->get(MW)),
         MTs = sqr(ival->get(MT)),
         MBs = sqr(ival->get(MB)),
         ALpi = ival->get(al)/Pi,
         ALS = ival->get(als);
  double cs = s/MTs;
  double rb = linex2d(axis1zhb, sizeof(axis1zhb)/sizeof(double),
                      axis2zhb, sizeof(axis2zhb)/sizeof(double),
		      s/MZs, MBs/MZs, &datarzhb[0][0]);
  double rt = linex1d(axis1zht, sizeof(axis1zht)/sizeof(double),
		      MZs/MTs, datar1zht) +
              linex1d(axis2zht, sizeof(axis2zht)/sizeof(double),
		      cs, datar2zht) - 4*Pi/3*(-Pi
		        *log(2*sqrt(abs((1-cs/4)/(1+cs/4)))) 
			+ 4*sqrt(fmax(1-cs/4,0)));
  /* the leading energy dependence near the ttbar threshold is provided via
     an analytical formula from NPB 347, 86, (18), to get a smoother dependence
     of the interpolation function on mt. */
  double ib = linex2d(axis1zhb, sizeof(axis1zhb)/sizeof(double),
                      axis2zhb, sizeof(axis2zhb)/sizeof(double),
		      s/MZs, MBs/MZs, &dataizhb[0][0]);
  double it = linex1d(axis2zht, sizeof(axis2zht)/sizeof(double),
		      cs, datai2zht) - ((cs>4) ? 
		       2*Pi/3*(Pi*Pi - 8*sqrt(cs/4-1)/(cs/4))/(cs/4) : 0);
  Cplx res = (20*MWs-11*MZs)*zh0(s/MZs)
               + (2*MWs-MZs/2)*(Cplx(rb,ib)+zh0(s/MZs))
	       + (8*MWs-5*MZs)*Cplx(rt,it);
  /* the dependence on the gauge couplings has been factored out */
  return(res * getalphas(s,MZs,ALS) 
  	  * ALpi/(sqrt(MWs*(MZs-MWs))*12.*Pi)); // exclude a factor e*Qf
}

// O(aas) corrections to the renormalized photon self-energy
Cplx aahataas(double s, const inval* ival)
{
  double MZs = sqr(ival->get(MZ)),
         MTs = sqr(ival->get(MT)),
         MBs = sqr(ival->get(MB)),
         ALpi = ival->get(al)/Pi,
         ALS = ival->get(als);
  double cs = s/MTs;
  double rb = linex2d(axis1aahb, sizeof(axis1aahb)/sizeof(double),
                      axis2aahb, sizeof(axis2aahb)/sizeof(double),
		      s/MZs, MBs/MZs, &dataraahb[0][0]);
  double rt = linex1d(axis1aaht, sizeof(axis1aaht)/sizeof(double),
		      MZs/MTs, datar1aaht) +
              linex1d(axis2aaht, sizeof(axis2aaht)/sizeof(double),
		      cs, datar2aaht) - 4*Pi/3*(-Pi
		        *log(2*sqrt(abs((1-cs/4)/(1+cs/4)))) 
			+ 4*sqrt(fmax(1-cs/4,0)));
  /* the leading energy dependence near the ttbar threshold is provided via
     an analytical formula from NPB 347, 86, (18), to get a smoother dependence
     of the interpolation function on mt. */
  double ib = linex2d(axis1aahb, sizeof(axis1aahb)/sizeof(double),
                      axis2aahb, sizeof(axis2aahb)/sizeof(double),
		      s/MZs, MBs/MZs, &dataiaahb[0][0]);
  double it = linex1d(axis2aaht, sizeof(axis2aaht)/sizeof(double),
		      cs, datai2aaht) - ((cs>4) ? 
		       2*Pi/3*(Pi*Pi - 8*sqrt(cs/4-1)/(cs/4))/(cs/4) : 0);
  Cplx res = 10*aah0(s/MZs) + (Cplx(rb,ib)+aah0(s/MZs)) + Cplx(rt,it);
  /* the dependence on the gauge couplings has been factored out */
  return(-res * getalphas(s,MZs,ALS) * ALpi*2/(3.*Pi) * s);
}

// O(aas) corrections to the renormalized Z-boson self-energy
Cplx zzhataas(double s, const inval* ival)
{
  double MZs = sqr(ival->get(MZ)),
         MWs = sqr(ival->get(MW)),
         MTs = sqr(ival->get(MT)),
         MBs = sqr(ival->get(MB)),
         ALpi = ival->get(al)/Pi,
         ALS = ival->get(als);
  if(abs(s-MZs) < 10)
  {
    Cplx z10 = zzhataas(MZs+10.01, ival);
    double ds = (s-MZs)/10.01;
    return(z10*ds*ds);
  }
  double cs = s/MTs;
  double rab= linex2d(axis1zahb, sizeof(axis1zahb)/sizeof(double),
                      axis2zahb, sizeof(axis2zahb)/sizeof(double),
		      s/MZs, MBs/MZs, &datarzahb[0][0]);
  double rvb= linex2d(axis1zvhb, sizeof(axis1zvhb)/sizeof(double),
                      axis2zvhb, sizeof(axis2zvhb)/sizeof(double),
		      s/MZs, MBs/MZs, &datarzvhb[0][0]);
  double rat= linex1d(axis1zaht, sizeof(axis1zaht)/sizeof(double),
		      MZs/MTs, datar0zaht)/cs +
              linex1d(axis1zaht, sizeof(axis1zaht)/sizeof(double),
		      MZs/MTs, datar1zaht) +
              linex1d(axis2zaht, sizeof(axis2zaht)/sizeof(double),
		      cs, datar2zaht) - 4*Pi/3*(-Pi
		        *log(2*sqrt(abs((1-cs/4)/(1+cs/4)))) 
			+ 4*sqrt(fmax(1-cs/4,0)));
  double rvt= linex1d(axis1zvht, sizeof(axis1zvht)/sizeof(double),
		      MZs/MTs, datar0zvht)/cs +
              linex1d(axis1zvht, sizeof(axis1zvht)/sizeof(double),
		      MZs/MTs, datar1zvht) +
              linex1d(axis2zvht, sizeof(axis2zvht)/sizeof(double),
		      cs, datar2zvht) - 4*Pi/3*(-Pi
		        *log(2*sqrt(abs((1-cs/4)/(1+cs/4)))) 
			+ 4*sqrt(fmax(1-cs/4,0)));
  /* the leading energy dependence near the ttbar threshold is provided via
     an analytical formula from NPB 347, 86, (18), to get a smoother dependence
     of the interpolation function on mt. */
  double iab= linex2d(axis1zahb, sizeof(axis1zahb)/sizeof(double),
                      axis2zahb, sizeof(axis2zahb)/sizeof(double),
		      s/MZs, MBs/MZs, &dataizahb[0][0]);
  double ivb= linex2d(axis1zvhb, sizeof(axis1zvhb)/sizeof(double),
                      axis2zvhb, sizeof(axis2zvhb)/sizeof(double),
		      s/MZs, MBs/MZs, &dataizvhb[0][0]);
  double iat= linex1d(axis2zaht, sizeof(axis2zaht)/sizeof(double),
		      cs, datai2zaht) - ((cs>4) ? 
		       2*Pi/3*(Pi*Pi - 8*sqrt(cs/4-1)/(cs/4))/(cs/4) : 0);
  double ivt= linex1d(axis2zvht, sizeof(axis2zvht)/sizeof(double),
		      cs, datai2zvht) - ((cs>4) ? 
		       2*Pi/3*(Pi*Pi - 8*sqrt(cs/4-1)/(cs/4))/(cs/4) : 0);
  Cplx res = -2*s*MWs*(10*MWs-11*MZs)*zvh0(s/MZs)
               - s*MWs*(2*MWs-MZs)*(Cplx(rvb,ivb)+zvh0(s/MZs))
	       - 2*s*MWs*(4*MWs-5*MZs)*Cplx(rvt,ivt)
             + (-44*s*zah0(s/MZs)
                - 5*s*(Cplx(rab,iab)+zah0(s/MZs))
	        - 25/8.*s*Cplx(rat,iat))
	       *MZs*MZs;
  /* the dependence on the gauge couplings has been factored out */
  return(res * getalphas(s,MZs,ALS) * ALpi/(12.*Pi) / MWs/(MZs-MWs));
}

Cplx mat_SMaas::coeffR(void) const
{
  if(iff > AXV || off > AXV)  // SCA and PSC contributions can only occur for 
    return 0;                 // Bhabha-like t-channel contributions, which does
                              // not have a s-channel resonance 

  double mz = ival->get(MZ),
	 gz = ival->get(GamZ);
  double QWe, QWf, IVe, IVf, rIAA, xI, IVe2, IVf2, rIAA2,
         zie0 = z0(it,iff,*ival), zjf0 = z0(ot,off,*ival),
	 iszp1 = isz1fp(*ival)+isz1bp(*ival), 
	 iszpp1 = isz1fpp(*ival)+isz1bpp(*ival);
  double MZs = mz*mz,
         MWs = sqr(ival->get(MW)),
	 ALpi = (ival->get(al))/Pi,
         ALS = ival->get(als);
#include "imaas.in"
  
  switch(iff)
  {
    case VEC: QWe = 1.-4*fabs(Qf[it])*realreg(SWi->result());
              IVe = (az0(it,*ival)*(iz1f(it,VEC,*ival)+iz1b(it,VEC,*ival)) 
	             - vz0(it,*ival)*(iz1f(it,AXV,*ival)+iz1b(it,AXV,*ival)))
		    /sqr(az0(it,*ival));
              IVe2 = (vg0(it,*ival)*izvaas)/az0(it,*ival);
	      break;
    case AXV: QWe = 1;  IVe = 0;  IVe2 = 0;
              break;
  }
  switch(off)
  {
    case VEC: QWf = 1.-4*fabs(Qf[ot])*realreg(SWo->result());
              IVf = (az0(ot,*ival)*(iz1f(ot,VEC,*ival)+iz1b(ot,VEC,*ival)) 
	             - vz0(ot,*ival)*(iz1f(ot,AXV,*ival)+iz1b(ot,AXV,*ival)))
		    /sqr(az0(ot,*ival));
              IVf2 = (vg0(ot,*ival)*izvaas)/az0(ot,*ival);
              break;
    case AXV: QWf = 1;  IVf = 0;  IVf2 = 0;
              break;
  }
  rIAA = (iz1f(it,AXV,*ival)+iz1b(it,AXV,*ival))/az0(it,*ival) 
  	 + (iz1f(ot,AXV,*ival)+iz1b(ot,AXV,*ival))/az0(ot,*ival) 
	 - iszp1;
  rIAA2 = -iszpaas;
  xI = (iz1fp(it,iff,*ival)+iz1bp(it,iff,*ival))/zie0 
       + (iz1fp(ot,off,*ival)+iz1bp(ot,off,*ival))/zjf0 - iszpp1/2;
  return(4*I3f[it]*I3f[ot]*sqrt(FAi->result()*FAo->result()) *
         (QWe*QWf*(Cplx(1 - rIAA*rIAA/2, rIAA + rIAA2) - iszp1*iszp1/2 
	           + bRaz1(it,ot,cost,*ival))
          + (QWe*IVf + QWf*IVe)*Cplx(-rIAA,1) - IVe*IVf
	  + (QWe*IVf2 + QWf*IVe2)*I)
	 + mz*gz*zie0*zjf0*xI);
}

Cplx mat_SMaas::resoffZaas(void) const
{
  double mz = ival->get(MZ);
  double gie0 = g0(it,iff,*ival), gjf0 = g0(ot,off,*ival),
         zie0 = z0(it,iff,*ival), zjf0 = z0(ot,off,*ival);
  Cplx zhat1 = zhataas(s,ival),
       aahat1 = aahataas(s,ival),
       zzhat1 = zzhataas(s,ival);
  if(it==ot)
    cerr << "The case of identical ini. and fin. flavors is not available in mat_SMaas yet!" << endl;  
  return((-(zie0*gjf0 + gie0*zjf0)*zhat1 - zie0*zjf0/(s-mz*mz)*zzhat1)/(s-mz*mz)
         - gie0*gjf0/(s*s)*aahat1);
}


Cplx mat_SMaas::result(void) const
{
  double mz = ival->get(MZ),
         gz = ival->get(GamZ);
  Cplx sminuss0(s - mz*mz, mz*gz);
  return(coeffR()/sminuss0 + resoffZ());
}

} // namespace

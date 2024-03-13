/*-----------------------------------------------------------------------------
xscnnlo.cc
Lisong Chen (lic114@pitt.edu), Ayres Freitas (afreitas@pitt.edu)
last revision: 14 Sep 2022
-------------------------------------------------------------------------------
matrix element class needed for description of cross-section at NNLO 
precision near Z pole and NLO away from Z pole
-----------------------------------------------------------------------------*/

#include "xscnnlo.h"
#include "ff0.h"
#include "ff.h"


Cplx mat_SMNNLO::coeffR(void) const
{
  double mz = ival->get(MZ),
	 gz = ival->get(GamZ);
  double QWe, QWf, IVe, IVf, rIAA, xI,
         zie0 = z0(it,iff,*ival), zjf0 = z0(ot,off,*ival),
	 iszp1 = isz1fp(*ival)+isz1bp(*ival), 
	 iszpp1 = isz1fpp(*ival)+isz1bpp(*ival);
  
  switch(iff)
  {
    case VEC: QWe = 1.-4*fabs(Qf[it])*realreg(SWi->result());
              IVe = (az0(it,*ival)*(iz1f(it,VEC,*ival)+iz1b(it,VEC,*ival)) 
	             - vz0(it,*ival)*(iz1f(it,AXV,*ival)+iz1b(it,AXV,*ival)))
		    /sqr(az0(it,*ival));
              break;
    case AXV: QWe = 1;  IVe = 0;
              break;
  }
  switch(off)
  {
    case VEC: QWf = 1.-4*fabs(Qf[ot])*realreg(SWo->result());
              IVf = (az0(ot,*ival)*(iz1f(ot,VEC,*ival)+iz1b(ot,VEC,*ival)) 
	             - vz0(ot,*ival)*(iz1f(ot,AXV,*ival)+iz1b(ot,AXV,*ival)))
		    /sqr(az0(ot,*ival));
              break;
    case AXV: QWf = 1;  IVf = 0;
              break;
  }
  rIAA = (iz1f(it,AXV,*ival)+iz1b(it,AXV,*ival))/az0(it,*ival) 
  	 + (iz1f(ot,AXV,*ival)+iz1b(ot,AXV,*ival))/az0(ot,*ival) 
	 - iszp1;
  xI = (iz1fp(it,iff,*ival)+iz1bp(it,iff,*ival))/zie0 
       + (iz1fp(ot,off,*ival)+iz1bp(ot,off,*ival))/zjf0 - iszpp1/2;
  return(4*I3f[it]*I3f[ot]*sqrt(FAi->result()*FAo->result()) *
         (QWe*QWf*(Cplx(1 - rIAA*rIAA/2, rIAA) - iszp1*iszp1/2 
	           /*+ bRaz1(it,ot,cost,*ival)*/)
// For use in POWHEG_EW: ^^^^^^^^^^^^^^^^^^^^ no gamZ box
          + (QWe*IVf + QWf*IVe)*Cplx(-rIAA,1) - IVe*IVf)
	 + mz*gz*zie0*zjf0*xI);
}

Cplx mat_SMNNLO::coeffS1f(void) const
{
  double mz = ival->get(MZ),
	 gz = ival->get(GamZ);
  double gie0 = g0(it,iff,*ival), gjf0 = g0(ot,off,*ival),
         zie0 = z0(it,iff,*ival), zjf0 = z0(ot,off,*ival);
  Cplx zpie1 = Cplx(rz1fp(it,iff,*ival), iz1fp(it,iff,*ival)),
       zpjf1 = Cplx(rz1fp(ot,off,*ival), iz1fp(ot,off,*ival)),
       gie1 = Cplx(rg1f(it,iff,*ival), ig1f(it,iff,*ival)),
       gjf1 = Cplx(rg1f(ot,off,*ival), ig1f(ot,off,*ival)),
       szpp1 = Cplx(rsz1fpp(*ival), isz1fpp(*ival)),
       sa1 = Cplx(rsg1f(*ival), isg1f(*ival)); 
  return(zie0*zpjf1 + zpie1*zjf0 - zie0*zjf0*szpp1/2
        + (gie0*gjf1 + gie1*gjf0 + gie0*gjf0*(/*I*gz/mz*/ - sa1/(mz*mz)))/(mz*mz));
// the gz/mz term is already in coeffS = S_SMLO   ^^^^^
}

Cplx mat_SMNNLO::coeffS1b(void) const
{
  double mz = ival->get(MZ),
	 gz = ival->get(GamZ);
  double gie0 = g0(it,iff,*ival), gjf0 = g0(ot,off,*ival),
         zie0 = z0(it,iff,*ival), zjf0 = z0(ot,off,*ival);
  Cplx zpie1 = Cplx(rz1bp(it,iff,*ival), iz1bp(it,iff,*ival)),
       zpjf1 = Cplx(rz1bp(ot,off,*ival), iz1bp(ot,off,*ival)),
       gie1 = Cplx(rg1b(it,iff,*ival), ig1b(it,iff,*ival)),
       gjf1 = Cplx(rg1b(ot,off,*ival), ig1b(ot,off,*ival)),
       szpp1 = Cplx(rsz1bpp(*ival), isz1bpp(*ival)),
       sa1 = Cplx(rsg1b(*ival), isg1b(*ival)); 
  return(zie0*zpjf1 + zpie1*zjf0 - zie0*zjf0*szpp1/2
        + (gie0*gjf1 + gie1*gjf0 + gie0*gjf0*(/*I*gz/mz*/ - sa1/(mz*mz)))/(mz*mz)
// the gz/mz term is already in coeffS = S_SMLO   ^^^^^
	+ B1(it,ot,iff,off,s,cost,*ival,0,0));
// Adjustment for use in POWHEG_EW:     ^^^ do not include gamgam and gamZ boxes
}

Cplx mat_SMNNLO::resoffZ1f(void) const
{
  double mz = ival->get(MZ);
  double gie0 = g0(it,iff,*ival), gjf0 = g0(ot,off,*ival),
         zie0 = z0(it,iff,*ival), zjf0 = z0(ot,off,*ival);
  Cplx zie1 = Cplx(rz1fs(it,iff,s,*ival), iz1fs(it,iff,s,*ival)),
       zie1z = Cplx(rz1f(it,iff,*ival), iz1f(it,iff,*ival)),
       zpie1z = Cplx(rz1fp(it,iff,*ival), iz1fp(it,iff,*ival)),
       zjf1 = Cplx(rz1fs(ot,off,s,*ival), iz1fs(ot,off,s,*ival)),
       zjf1z = Cplx(rz1f(ot,off,*ival), iz1f(ot,off,*ival)),
       zpjf1z = Cplx(rz1fp(ot,off,*ival), iz1fp(ot,off,*ival)),
       gie1 = Cplx(rg1fs(it,iff,s,*ival), ig1fs(it,iff,s,*ival)),
       gie1z = Cplx(rg1f(it,iff,*ival), ig1f(it,iff,*ival)),
       gjf1 = Cplx(rg1fs(ot,off,s,*ival), ig1fs(ot,off,s,*ival)),
       gjf1z = Cplx(rg1f(ot,off,*ival), ig1f(ot,off,*ival)),
       sz1 = Cplx(rsz1fs(s,*ival), isz1fs(s,*ival)),
       sz1z = Cplx(rsz1f(*ival), isz1f(*ival)),
       szp1z = Cplx(rsz1fp(*ival), isz1fp(*ival)),
       szpp1z = Cplx(rsz1fpp(*ival), isz1fpp(*ival)),
       sa1 = Cplx(rsg1fs(s,*ival), isg1fs(s,*ival)), 
       sa1z = Cplx(rsg1f(*ival), isg1f(*ival)); 
  Cplx Rp = -zie0*zjf0*sz1z,
       R = zie0*zjf1z + zie1z*zjf0 - zie0*zjf0*szp1z,
       S = (zie0*zpjf1z + zpie1z*zjf0 - zie0*zjf0*szpp1z/2
          + (gie0*gjf1z + gie1z*gjf0 - gie0*gjf0*sa1z/(mz*mz))/(mz*mz)),
       mats1 = ((zie0*zjf1 + zie1*zjf0 - zie0*zjf0*sz1/(s-mz*mz))/(s-mz*mz)
          + (gie0*gjf1 + gie1*gjf0 - gie0*gjf0*sa1/s)/s);
  return(mats1 - ((R + Rp/(s-mz*mz))/(s-mz*mz) + S));
}

Cplx mat_SMNNLO::resoffZ1b(void) const
{
  double mz = ival->get(MZ);
  double gie0 = g0(it,iff,*ival), gjf0 = g0(ot,off,*ival),
         zie0 = z0(it,iff,*ival), zjf0 = z0(ot,off,*ival);
  Cplx zie1 = Cplx(rz1bs(it,iff,s,*ival), iz1bs(it,iff,s,*ival)),
       zie1z = Cplx(rz1b(it,iff,*ival), iz1b(it,iff,*ival)),
       zpie1z = Cplx(rz1bp(it,iff,*ival), iz1bp(it,iff,*ival)),
       zjf1 = Cplx(rz1bs(ot,off,s,*ival), iz1bs(ot,off,s,*ival)),
       zjf1z = Cplx(rz1b(ot,off,*ival), iz1b(ot,off,*ival)),
       zpjf1z = Cplx(rz1bp(ot,off,*ival), iz1bp(ot,off,*ival)),
       gie1 = Cplx(rg1bs(it,iff,s,*ival), ig1bs(it,iff,s,*ival)),
       gie1z = Cplx(rg1b(it,iff,*ival), ig1b(it,iff,*ival)),
       gjf1 = Cplx(rg1bs(ot,off,s,*ival), ig1bs(ot,off,s,*ival)),
       gjf1z = Cplx(rg1b(ot,off,*ival), ig1b(ot,off,*ival)),
       sz1 = Cplx(rsz1bs(s,*ival), isz1bs(s,*ival)),
       sz1z = Cplx(rsz1b(*ival), isz1b(*ival)),
       szp1z = Cplx(rsz1bp(*ival), isz1bp(*ival)),
       szpp1z = Cplx(rsz1bpp(*ival), isz1bpp(*ival)),
       sa1 = Cplx(rsg1bs(s,*ival), isg1bs(s,*ival)), 
       sa1z = Cplx(rsg1b(*ival), isg1b(*ival)); 
  Cplx Rp = -zie0*zjf0*sz1z,
       R = zie0*zjf1z + zie1z*zjf0 + zie0*zjf0*(-szp1z /*+ bRaz1(it,ot,cost,*ival)*/),
// For use in POWHEG_EW: no gamZ box                     ^^^^^^^^^^^^^^^^^^^^^^^^^ 
       S = (zie0*zpjf1z + zpie1z*zjf0 - zie0*zjf0*szpp1z/2
          + (gie0*gjf1z + gie1z*gjf0 - gie0*gjf0*sa1z/(mz*mz))/(mz*mz)
	  + B1(it,ot,iff,off,s,cost,*ival,0,0,1e-12)),
// Adjustment for use in POWHEG_EW:       ^^^ do not include gamgam and gamZ boxes
       mats1 = ((zie0*zjf1 + zie1*zjf0 - zie0*zjf0*sz1/(s-mz*mz))/(s-mz*mz)
          + (gie0*gjf1 + gie1*gjf0 - gie0*gjf0*sa1/s)/s
	  + B1s(it,ot,iff,off,s,cost,*ival,0,0));
// Adjustment for use in POWHEG_EW:        ^^^ do not include gamgam and gamZ boxes
  return(mats1 - ((R + Rp/(s-mz*mz))/(s-mz*mz) + S));
}

Cplx mat_SMNNLO::result(void) const
{
  double mz = ival->get(MZ),
         gz = ival->get(GamZ);
  Cplx sminuss0(s - mz*mz, mz*gz);
  return(coeffR()/sminuss0 + coeffS() + coeffSp()*sminuss0 + resoffZ());
}

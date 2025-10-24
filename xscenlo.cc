/*-----------------------------------------------------------------------------
xscenlo.cc
Ayres Freitas (afreitas@pitt.edu)
last revision: 24 Oct 2025
-------------------------------------------------------------------------------
matrix element at EXACTLY NLO (no products of 1-loop terms), for testing
-----------------------------------------------------------------------------*/

#include "xscenlo.h"
//#include "ff0.h"
#include "ff.h"

namespace griffin {

Cplx mat_SMeNLO::coeffR(void) const
{
  double mz = ival->get(MZ),
	 gz = ival->get(GamZ),
	 mw = ival->get(MW);
  double sw0 = 1-sqr(mw/mz);
  double QWe, QWf, QWe0, QWf0, IVe, IVf, rIAA, xI,
         zie0 = z0(it,iff,*ival), zjf0 = z0(ot,off,*ival),
	 iszp1 = isz1fp(*ival)+isz1bp(*ival), 
	 iszpp1 = isz1fpp(*ival)+isz1bpp(*ival);
  
  switch(iff)
  {
    case VEC: QWe = 1.-4*fabs(Qf[it])*realreg(SWi->result());
              QWe0 = 1.-4*fabs(Qf[it])*sw0;
              IVe = (az0(it,*ival)*(iz1f(it,VEC,*ival)+iz1b(it,VEC,*ival)) 
	             - vz0(it,*ival)*(iz1f(it,AXV,*ival)+iz1b(it,AXV,*ival)))
		    /sqr(az0(it,*ival));
              break;
    case AXV: QWe = QWe0 = 1;  IVe = 0;
              break;
  }
  switch(off)
  {
    case VEC: QWf = 1.-4*fabs(Qf[ot])*realreg(SWo->result());
              QWf0 = 1.-4*fabs(Qf[ot])*sw0;
              IVf = (az0(ot,*ival)*(iz1f(ot,VEC,*ival)+iz1b(ot,VEC,*ival)) 
	             - vz0(ot,*ival)*(iz1f(ot,AXV,*ival)+iz1b(ot,AXV,*ival)))
		    /sqr(az0(ot,*ival));
              break;
    case AXV: QWf = QWf0 = 1;  IVf = 0;
              break;
  }
  rIAA = (iz1f(it,AXV,*ival)+iz1b(it,AXV,*ival))/az0(it,*ival) 
  	 + (iz1f(ot,AXV,*ival)+iz1b(ot,AXV,*ival))/az0(ot,*ival) 
	 - iszp1;
  xI = (iz1fp(it,iff,*ival)+iz1bp(it,iff,*ival))/zie0 
       + (iz1fp(ot,off,*ival)+iz1bp(ot,off,*ival))/zjf0 - iszpp1/2;
//  return(4*I3f[it]*I3f[ot]*sqrt(FAi->result()*FAo->result()) *
//         (QWe*QWf*(Cplx(1 - rIAA*rIAA/2, rIAA) - iszp1*iszp1/2 
//	           /*+ bRaz1(it,ot,cost,*ival)*/)
//          + (QWe*IVf + QWf*IVe)*Cplx(-rIAA,1) - IVe*IVf)
//	 + mz*gz*zie0*zjf0*xI);
  double FAi0 = sqr(az0(it,*ival)), FAo0 = sqr(az0(ot,*ival));
  Cplx delFAi = FAi->result() - FAi0,
       delFAo = FAo->result() - FAo0;
  double delQWe = QWe-QWe0,
	 delQWf = QWf-QWf0;
  return(4*I3f[it]*I3f[ot]*sqrt(FAi0*FAo0)*(QWe0*QWf0 *
       (Cplx(1, rIAA) + (delFAi/FAi0 + delFAo/FAo0)/2. 
         + delQWe/QWe0 + delQWf/QWf0) + I*(QWe0*IVf + QWf0*IVe))
	 /*+ mz*gz*zie0*zjf0*xI*/);
}

Cplx mat_SMeNLO::resoffZ1f(void) const
{
  double mz = ival->get(MZ), stmp = s;
  double gie0 = g0(it,iff,*ival), gjf0 = g0(ot,off,*ival),
         zie0 = z0(it,iff,*ival), zjf0 = z0(ot,off,*ival);
  if(abs(1-s/(mz*mz)) < 1e-5)
    stmp *= 1+1e-8;
  Cplx zie1 = Cplx(rz1fs(it,iff,stmp,*ival), iz1fs(it,iff,stmp,*ival)),
       zie1z = Cplx(rz1f(it,iff,*ival), iz1f(it,iff,*ival)),
       zjf1 = Cplx(rz1fs(ot,off,stmp,*ival), iz1fs(ot,off,stmp,*ival)),
       zjf1z = Cplx(rz1f(ot,off,*ival), iz1f(ot,off,*ival)),
       gie1 = Cplx(rg1fs(it,iff,s,*ival), ig1fs(it,iff,s,*ival)),
       gjf1 = Cplx(rg1fs(ot,off,s,*ival), ig1fs(ot,off,s,*ival)),
       sz1 = Cplx(rsz1fs(s,*ival), isz1fs(s,*ival)),
       sz1z = Cplx(rsz1f(*ival), isz1f(*ival)),
       szp1z = Cplx(rsz1fp(*ival), isz1fp(*ival)),
       sa1 = Cplx(rsg1fs(s,*ival), isg1fs(s,*ival)); 
  if(abs(1-s/(mz*mz)) < 1e-5)
    sz1 = sz1z = szp1z = 0;
  Cplx Rp = -zie0*zjf0*sz1z,
       R = zie0*zjf1z + zie1z*zjf0 - zie0*zjf0*szp1z,
       mats1 = ((zie0*zjf1 + zie1*zjf0 - zie0*zjf0*sz1/(stmp-mz*mz))/(stmp-mz*mz)
          + (gie0*gjf1 + gie1*gjf0 - gie0*gjf0*sa1/s)/s);
  return(mats1 - ((R + Rp/(stmp-mz*mz))/(stmp-mz*mz) /*+ S*/));
}

Cplx mat_SMeNLO::resoffZ1b(void) const
{
  double mz = ival->get(MZ), stmp = s;
  double gie0 = g0(it,iff,*ival), gjf0 = g0(ot,off,*ival),
         zie0 = z0(it,iff,*ival), zjf0 = z0(ot,off,*ival);
  if(abs(1-s/(mz*mz)) < 1e-5)
    stmp *= 1+1e-8;
  Cplx zie1 = Cplx(rz1bs(it,iff,stmp,*ival), iz1bs(it,iff,stmp,*ival)),
       zie1z = Cplx(rz1b(it,iff,*ival), iz1b(it,iff,*ival)),
       zjf1 = Cplx(rz1bs(ot,off,stmp,*ival), iz1bs(ot,off,stmp,*ival)),
       zjf1z = Cplx(rz1b(ot,off,*ival), iz1b(ot,off,*ival)),
       gie1 = Cplx(rg1bs(it,iff,s,*ival), ig1bs(it,iff,s,*ival)),
       gjf1 = Cplx(rg1bs(ot,off,s,*ival), ig1bs(ot,off,s,*ival)),
       sz1 = Cplx(rsz1bs(s,*ival), isz1bs(s,*ival)),
       sz1z = Cplx(rsz1b(*ival), isz1b(*ival)),
       szp1z = Cplx(rsz1bp(*ival), isz1bp(*ival)),
       sa1 = Cplx(rsg1bs(s,*ival), isg1bs(s,*ival)); 
  if(abs(1-s/(mz*mz)) < 1e-5)
    sz1 = sz1z = szp1z = 0;
  Cplx Rp = -zie0*zjf0*sz1z,
       R = zie0*zjf1z + zie1z*zjf0 + zie0*zjf0*(-szp1z /*+ bRaz1(it,ot,cost,*ival)*/),
// For use in POWHEG_EW: no gamZ box                     ^^^^^^^^^^^^^^^^^^^^^^^^^ 
       mats1 = ((zie0*zjf1 + zie1*zjf0 - zie0*zjf0*sz1/(stmp-mz*mz))/(stmp-mz*mz)
          + (gie0*gjf1 + gie1*gjf0 - gie0*gjf0*sa1/s)/s
	  + B1s(it,ot,iff,off,stmp,-stmp/2*(1-cost),*ival,0,0,1));
// Adjustment for POWHEG_EW: no gamgam/gamZ boxes         ^^^ 
  return(mats1 - ((R + Rp/(stmp-mz*mz))/(stmp-mz*mz) /*+ S*/));
}

Cplx mat_SMeNLO::result(void) const
{
  double mz = ival->get(MZ),
         gz = ival->get(GamZ);
  Cplx sminuss0(s - mz*mz, mz*gz);
  return(coeffR()/sminuss0 /*+ coeffS() + coeffSp()*sminuss0*/ + resoffZ());
}

} // namespace

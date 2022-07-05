/*-----------------------------------------------------------------------------
xscnnlo.cc
Lisong Chen (lic114@pitt.edu), Ayres Freitas (afreitas@pitt.edu)
last revision: 16 Dec 2021
-------------------------------------------------------------------------------
matrix element classes needed for description of cross-section at NNLO 
precision near Z pole
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
	           + bRaz1(it,ot,cost,*ival))
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
	+ B1(it,ot,iff,off,mz*mz,cost,*ival,1,1));
}

Cplx mat_SMNNLO::result(void) const
{
  double mz = ival->get(MZ),
         gz = ival->get(GamZ);
  Cplx sminuss0(s - mz*mz, mz*gz);
  return(coeffR()/sminuss0 + coeffS() + coeffSp()*sminuss0);
}


Cplx msq_SMNNLO::result(void) const
{
  double mz = ival->get(MZ),
         gz = ival->get(GamZ);
  double gie0 = g0(it,if1,*ival), gke0 = g0(it,if2,*ival),
         gjf0 = g0(ot,of1,*ival), glf0 = g0(ot,of2,*ival),
	 zie0 = z0(it,if1,*ival), zke0 = z0(it,if2,*ival),
         zjf0 = z0(ot,of1,*ival), zlf0 = z0(ot,of2,*ival);
  double zpie1 = rz1fp(it,if1,*ival)+rz1bp(it,if1,*ival),
	 zpke1 = rz1fp(it,if2,*ival)+rz1bp(it,if2,*ival),
	 zpjf1 = rz1fp(ot,of1,*ival)+rz1bp(ot,of1,*ival),
	 zplf1 = rz1fp(ot,of2,*ival)+rz1bp(ot,of2,*ival),
	 iszp1 = isz1fp(*ival)+isz1bp(*ival),
	 szpp1 = rsz1fpp(*ival)+rsz1bpp(*ival);
  Cplx zie1 = Cplx(rz1f(it,if1,*ival)+rz1b(it,if1,*ival), 
                   iz1f(it,if1,*ival)+iz1b(it,if1,*ival)),
       zke1 = Cplx(rz1f(it,if2,*ival)+rz1b(it,if2,*ival), 
                   iz1f(it,if2,*ival)+iz1b(it,if2,*ival)),
       zjf1 = Cplx(rz1f(ot,of1,*ival)+rz1b(ot,of1,*ival), 
                   iz1f(ot,of1,*ival)+iz1b(ot,of1,*ival)),
       zlf1 = Cplx(rz1f(ot,of2,*ival)+rz1b(ot,of2,*ival), 
                   iz1f(ot,of2,*ival)+iz1b(ot,of2,*ival)),
       gie1 = Cplx(rg1f(it,if1,*ival)+rg1b(it,if1,*ival), 
                   ig1f(it,if1,*ival)+ig1b(it,if1,*ival)),
       gke1 = Cplx(rg1f(it,if2,*ival)+rg1b(it,if2,*ival), 
                   ig1f(it,if2,*ival)+ig1b(it,if2,*ival)),
       gjf1 = Cplx(rg1f(ot,of1,*ival)+rg1b(ot,of1,*ival), 
                   ig1f(ot,of1,*ival)+ig1b(ot,of1,*ival)),
       glf1 = Cplx(rg1f(ot,of2,*ival)+rg1b(ot,of2,*ival), 
                   ig1f(ot,of2,*ival)+ig1b(ot,of2,*ival)),
       szp1 = Cplx(rsz1fp(*ival)+rsz1bp(*ival), iszp1),
       sa1 = Cplx(rsg1f(*ival)+rsg1b(*ival), isg1f(*ival)+isg1b(*ival)); 
  double QWie, QWjf, QWke, QWlf, IVie, IVjf, IVke, IVlf,
         Uijkl, Vijkl, Xijkl;
  Cplx Yijkl; 

  switch(if1)
  {
    case VEC: QWie = 1.-4*fabs(Qf[it])*realreg(SWi->result());
              IVie = (az0(it,*ival)*(iz1f(it,VEC,*ival)+iz1b(it,VEC,*ival)) 
	              - vz0(it,*ival)*(iz1f(it,AXV,*ival)+iz1b(it,AXV,*ival)))
		     /sqr(az0(it,*ival));
              break;
    case AXV: QWie = 1;  IVie = 0;
              break;
  }
  switch(if2)
  {
    case VEC: QWke = 1.-4*fabs(Qf[it])*realreg(SWi->result());
              IVke = (az0(it,*ival)*(iz1f(it,VEC,*ival)+iz1b(it,VEC,*ival)) 
	              - vz0(it,*ival)*(iz1f(it,AXV,*ival)+iz1b(it,AXV,*ival)))
		     /sqr(az0(it,*ival));
              break;
    case AXV: QWke = 1;  IVke = 0;
              break;
  }
  switch(of1)
  {
    case VEC: QWjf = 1.-4*fabs(Qf[ot])*realreg(SWo->result());
              IVjf = (az0(ot,*ival)*(iz1f(ot,VEC,*ival)+iz1b(ot,VEC,*ival)) 
	              - vz0(ot,*ival)*(iz1f(ot,AXV,*ival)+iz1b(ot,AXV,*ival)))
		     /sqr(az0(ot,*ival));
              break;
    case AXV: QWjf = 1;  IVjf = 0;
              break;
  }
  switch(of2)
  {
    case VEC: QWlf = 1.-4*fabs(Qf[ot])*realreg(SWo->result());
              IVlf = (az0(ot,*ival)*(iz1f(ot,VEC,*ival)+iz1b(ot,VEC,*ival)) 
	              - vz0(ot,*ival)*(iz1f(ot,AXV,*ival)+iz1b(ot,AXV,*ival)))
		     /sqr(az0(ot,*ival));
              break;
    case AXV: QWlf = 1;  IVlf = 0;
              break;
  }
  Uijkl = QWie*QWke*IVjf*IVlf + QWjf*QWlf*IVie*IVke
          + (QWie*IVke - QWke*IVie)*(QWlf*IVjf - QWjf*IVlf);
  Xijkl = zie0*zjf0*zke0*zplf1 + zie0*zjf0*zpke1*zlf0
          + zie0*zpjf1*zke0*zlf0 + zpie1*zjf0*zke0*zlf0
	  - zie0*zjf0*zke0*zlf0*szpp1;
  Yijkl = zie0*zjf0*(gke0*glf1 + gke1*glf0) 
          + gke0*glf0*conj(zie0*zjf1 + zie1*zjf0)
	  + zke0*zlf0*(gie0*gjf1 + gie1*gjf0) 
          + gie0*gjf0*conj(zke0*zlf1 + zke1*zlf0)
	  + (zie0*zjf0*gke0*glf0 + zke0*zlf0*gie0*gjf0)
	    *(1 - conj(szp1) - sa1/(mz*mz))
	  + mz*mz*(zie0*zjf0*B1(it,ot,if2,of2,mz*mz,cost,*ival,1,1)
	         + zke0*zlf0*B1(it,ot,if1,of1,mz*mz,cost,*ival,1,1));
  Vijkl = gie0*gjf0*(gke0*glf0/2 - zke0*zlf0)
          + gke0*glf0*(gie0*gjf0/2 - zie0*zjf0);
  return( (FAi->result()*FAo->result()*(QWie*QWjf*QWke*QWlf
  	     *(1 - iszp1*iszp1 + 2*bRaz1(it,ot,cost,*ival)) + Uijkl)
	   - gz/mz*imag(Yijkl) + gz*gz/(mz*mz)*gie0*gjf0*gke0*glf0
	   + (s - mz*mz) * (Xijkl + real(Yijkl)/(mz*mz))
	   + sqr(s/(mz*mz) - 1) * Vijkl) / (sqr(s - mz*mz) + sqr(mz*gz)) );
}

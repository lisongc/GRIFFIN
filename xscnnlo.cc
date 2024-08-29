/*-----------------------------------------------------------------------------
xscnnlo.cc
Lisong Chen (lic114@pitt.edu), Ayres Freitas (afreitas@pitt.edu)
last revision: 15 Aug 2024
-------------------------------------------------------------------------------
matrix element class needed for description of cross-section at NNLO 
precision near Z pole and NLO away from Z pole
-----------------------------------------------------------------------------*/

#include "xscnnlo.h"
#include "ff0.h"
#include "ff.h"

namespace griffin {

Cplx mat_SMNNLO::coeffR(void) const
{
  if(iff > AXV || off > AXV)  // SCA and PSC contributions can only occur for 
    return 0;                 // Bhabha-like t-channel contributions, which does
                              // not have a s-channel resonance 
  
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
	+ B1(it,ot,iff,off,s,cost,*ival,1,1));
}

Cplx mat_SMNNLO::resoffZ1f(void) const
{
  double mz = ival->get(MZ);
  double gie0 = g0(it,iff,*ival), gjf0 = g0(ot,off,*ival),
         zie0 = z0(it,iff,*ival), zjf0 = z0(ot,off,*ival);
  Cplx zie1 = Cplx(rz1fs(it,iff,s,*ival), iz1fs(it,iff,s,*ival)),
       zie1z = Cplx(rz1f(it,iff,*ival), iz1f(it,iff,*ival)),
       zjf1 = Cplx(rz1fs(ot,off,s,*ival), iz1fs(ot,off,s,*ival)),
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
       mats1 = ((zie0*zjf1 + zie1*zjf0 - zie0*zjf0*sz1/(s-mz*mz))/(s-mz*mz)
          + (gie0*gjf1 + gie1*gjf0 - gie0*gjf0*sa1/s)/s);
  if(it==ot)
  {
    double t = -s/2*(1-cost);
    double g0v = g0(it,VEC,*ival), g0a = g0(it,AXV,*ival),
    	   z0v = z0(it,VEC,*ival), z0a = z0(it,AXV,*ival);
    Cplx z1v = Cplx(rz1fs(it,VEC,t,*ival), iz1fs(it,VEC,t,*ival)),
         z1a = Cplx(rz1fs(it,AXV,t,*ival), iz1fs(it,AXV,t,*ival)),
	 g1v = Cplx(rg1fs(it,VEC,t,*ival), ig1fs(it,VEC,t,*ival)),
         g1a = Cplx(rg1fs(it,AXV,t,*ival), ig1fs(it,AXV,t,*ival)),
	 sz1 = Cplx(rsz1fs(t,*ival), isz1fs(t,*ival)),
	 sa1 = Cplx(rsg1fs(t,*ival), isg1fs(t,*ival));
    Cplx mats2;
    if(iff==off)
    {
      if(iff==VEC || iff==AXV)
        mats2 = ((2*(z0v*z1v+z0a*z1a) - (z0v*z0v+z0a*z0a)*sz1/(t-mz*mz))/(t-mz*mz)
                + (2*(g0v*g1v+g0a*g1a) - (g0v*g0v+g0a*g0a)*sa1/t)/t)/2;
      else
        mats2 = ((2*(z0v*z1v-z0a*z1a) - (z0v*z0v-z0a*z0a)*sz1/(t-mz*mz))/(t-mz*mz)
                + (2*(g0v*g1v-g0a*g1a) - (g0v*g0v-g0a*g0a)*sa1/t)/t);
    }
    else
    {
      if(iff==VEC || iff==AXV)
        mats2 = ((2*(z0v*z1a+z1v*z0a) - 2*z0v*z0a*sz1/(t-mz*mz))/(t-mz*mz)
                + (2*(g0v*g1a+g1v*g0a) - 2*g0v*g0a*sa1/t)/t)/2;
      else
        mats2 = 0;
    }
    mats1 += mats2;
  }
  return(mats1 - ((R + Rp/(s-mz*mz))/(s-mz*mz)));
}

Cplx mat_SMNNLO::resoffZ1b(void) const
{
  double mz = ival->get(MZ);
  double gie0 = g0(it,iff,*ival), gjf0 = g0(ot,off,*ival),
         zie0 = z0(it,iff,*ival), zjf0 = z0(ot,off,*ival);
  Cplx zie1 = Cplx(rz1bs(it,iff,s,*ival), iz1bs(it,iff,s,*ival)),
       zie1z = Cplx(rz1b(it,iff,*ival), iz1b(it,iff,*ival)),
       zjf1 = Cplx(rz1bs(ot,off,s,*ival), iz1bs(ot,off,s,*ival)),
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
       R = zie0*zjf1z + zie1z*zjf0 + zie0*zjf0*(-szp1z + bRaz1(it,ot,cost,*ival)),
       mats1 = ((zie0*zjf1 + zie1*zjf0 - zie0*zjf0*sz1/(s-mz*mz))/(s-mz*mz)
          + (gie0*gjf1 + gie1*gjf0 - gie0*gjf0*sa1/s)/s
	  + B1s(it,ot,iff,off,s,-s/2*(1-cost),*ival,1,1,1))
	  + B1s0(it,ot,iff,off,s,cost,*ival);
  if(it==ot)
  {
    double t = -s/2*(1-cost);
    double g0v = g0(it,VEC,*ival), g0a = g0(it,AXV,*ival),
    	   z0v = z0(it,VEC,*ival), z0a = z0(it,AXV,*ival);
    Cplx z1v = Cplx(rz1bs(it,VEC,t,*ival), iz1bs(it,VEC,t,*ival)),
         z1a = Cplx(rz1bs(it,AXV,t,*ival), iz1bs(it,AXV,t,*ival)),
	 g1v = Cplx(rg1bs(it,VEC,t,*ival), ig1bs(it,VEC,t,*ival)),
         g1a = Cplx(rg1bs(it,AXV,t,*ival), ig1bs(it,AXV,t,*ival)),
	 sz1 = Cplx(rsz1bs(t,*ival), isz1bs(t,*ival)),
	 sa1 = Cplx(rsg1bs(t,*ival), isg1bs(t,*ival));
    Cplx mats2, box2;
    if(iff==off)
    {
      if(iff==VEC || iff==AXV)  // iff=off=VEC or AXV
      {
        mats2 = ((2*(z0v*z1v+z0a*z1a) - (z0v*z0v+z0a*z0a)*sz1/(t-mz*mz))/(t-mz*mz)
                + (2*(g0v*g1v+g0a*g1a) - (g0v*g0v+g0a*g0a)*sa1/t)/t)/2;
        box2 = (B1s(it,ot,VEC,VEC,t,s,*ival,1,1,0)
	       +B1s(it,ot,AXV,AXV,t,s,*ival,1,1,0))/2;
      }
      else  // iff=off=SCA or PSC
      {
        mats2 = ((2*(z0v*z1v-z0a*z1a) - (z0v*z0v-z0a*z0a)*sz1/(t-mz*mz))/(t-mz*mz)
                + (2*(g0v*g1v-g0a*g1a) - (g0v*g0v-g0a*g0a)*sa1/t)/t);
        box2 = (B1s(it,ot,VEC,VEC,t,s,*ival,1,1,0)
	       -B1s(it,ot,AXV,AXV,t,s,*ival,1,1,0));
      }
    }
    else
    {
      if(iff==VEC || iff==AXV)  // iff=VEC and off=AXV, or vice versa
      {
        mats2 = ((2*(z0v*z1a+z1v*z0a) - 2*z0v*z0a*sz1/(t-mz*mz))/(t-mz*mz)
                + (2*(g0v*g1a+g1v*g0a) - 2*g0v*g0a*sa1/t)/t)/2;
        box2 = (B1s(it,ot,AXV,VEC,t,s,*ival,1,1,0)
	       +B1s(it,ot,VEC,AXV,t,s,*ival,1,1,0))/2;
      }
      else  // iff=SCA and off=PSC, or vice versa
      {
        mats2 = 0;
        box2 = (B1s(it,ot,AXV,VEC,t,s,*ival,1,1,0)
	       -B1s(it,ot,VEC,AXV,t,s,*ival,1,1,0));
	       // apparently this is always 0 for massless external fermions, 
	       // because some chirality flip is always required	       
      }
    }
    mats1 += mats2 + box2;
  }
  return(mats1 - ((R + Rp/(s-mz*mz))/(s-mz*mz)));
}

Cplx mat_SMNNLO::result(void) const
{
  double mz = ival->get(MZ),
         gz = ival->get(GamZ);
  Cplx sminuss0(s - mz*mz, mz*gz);
  return(coeffR()/sminuss0 + resoffZ());
}

} // namespace

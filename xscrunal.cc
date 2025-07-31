/*-----------------------------------------------------------------------------
xscrunal.cc
Ayres Freitas (afreitas@pitt.edu)
last revision: 30 July 2025
-------------------------------------------------------------------------------
similar to mat_SMNNLO, but with running alpha(s/t) in the s/t-channel QED 
amplitudes; alpha(s/t) are supplied as external parameters
-----------------------------------------------------------------------------*/

#include "xscrunal.h"
#include "ff0.h"
#include "ff.h"

namespace griffin {

double g0s(int type, int formt, const inval& input)
{
  if(formt==VEC)
  { return -sqrt(4*Pi * input.get(alQs))*Qf[type]; }
  else
  { return 0; }
}

double g0t(int type, int formt, const inval& input)
{
  if(formt==VEC)
  { return -sqrt(4*Pi * input.get(alQt))*Qf[type]; }
  else
  { return 0; }
}


Cplx mat_SMNNLOrunal::resoffZ0(void) const
{
  double mz = ival->get(MZ);
  double gi = g0s(it,iff,*ival), gf = g0s(ot,off,*ival);
  double res = gi*gf/s;
  if(it==ot)
  {
    double t = -s/2*(1-cost);
    double g0v = g0t(it,VEC,*ival), g0a = g0t(it,AXV,*ival),
    	   z0v = z0(it,VEC,*ival), z0a = z0(it,AXV,*ival);
    if(iff==off)
    {
      if(iff==VEC || iff==AXV)  // iff=off=VEC or AXV
        res += ((z0v*z0v+z0a*z0a)/(t-mz*mz) + (g0v*g0v+g0a*g0a)/t)/2;
      else                      // iff=off=SCA or PSC
        res += ((z0v*z0v-z0a*z0a)/(t-mz*mz) + (g0v*g0v-g0a*g0a)/t);
    }
    else
    {
      if(iff==VEC || iff==AXV)  // iff=VEC and off=AXV, or vice versa
        res += (2*z0v*z0a/(t-mz*mz) + 2*g0v*g0a/t)/2;
    }
  }
  return(res);
}

Cplx mat_SMNNLOrunal::resoffZ1f(void) const
{
  double mz = ival->get(MZ);
  double gie0 = g0(it,iff,*ival), gjf0 = g0(ot,off,*ival),
         zie0 = z0(it,iff,*ival), zjf0 = z0(ot,off,*ival);
  Cplx zie1 = Cplx(rz1fs(it,iff,s,*ival), iz1fs(it,iff,s,*ival)),
       zie1z = Cplx(rz1f(it,iff,*ival), iz1f(it,iff,*ival)),
       zjf1 = Cplx(rz1fs(ot,off,s,*ival), iz1fs(ot,off,s,*ival)),
       zjf1z = Cplx(rz1f(ot,off,*ival), iz1f(ot,off,*ival)),
       gie1 = Cplx(0, ig1fs(it,iff,s,*ival)),
       gjf1 = Cplx(0, ig1fs(ot,off,s,*ival)),
       sz1 = Cplx(rsz1fs(s,*ival), isz1fs(s,*ival)),
       sz1z = Cplx(rsz1f(*ival), isz1f(*ival)),
       szp1z = Cplx(rsz1fp(*ival), isz1fp(*ival)),
       sa1 = Cplx(0, isg1fs(s,*ival)); 
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
	 g1v = Cplx(0, ig1fs(it,VEC,t,*ival)),
         g1a = Cplx(0, ig1fs(it,AXV,t,*ival)),
	 sz1 = Cplx(rsz1fs(t,*ival), isz1fs(t,*ival)),
	 sa1 = Cplx(0, isg1fs(t,*ival));
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

Cplx mat_SMNNLOrunal::result(void) const
{
  double mz = ival->get(MZ),
         gz = ival->get(GamZ);
  Cplx sminuss0(s - mz*mz, mz*gz);
  return(coeffR()/sminuss0 + resoffZ());
}

} // namespace

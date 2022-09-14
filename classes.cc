/*-----------------------------------------------------------------------------
classes.cc
Lisong Chen (lic114@pitt.edu), Ayres Freitas (afreitas@pitt.edu)
last revision: 10 Feb 2022
-------------------------------------------------------------------------------
basic classes for form factors and matrix elements, including SM LO predictions
-----------------------------------------------------------------------------*/

#include "classes.h"
#include "ff0.h"

// the base classes compute results at LO:

Cplx matel::coeffR(void) const
{
  Cplx zi, zf;
  
  switch(iff)
  {
    case VEC: zi = 2*I3f[it]*sqrt(FAi->result())
    			*(1.-4*fabs(Qf[it])*realreg(SWi->result()));
              break;
    case AXV: zi = 2*I3f[it]*sqrt(FAi->result());
              break;
  }
  switch(off)
  {
    case VEC: zf = 2*I3f[ot]*sqrt(FAo->result())
    			*(1.-4*fabs(Qf[ot])*realreg(SWo->result()));
              break;
    case AXV: zf = 2*I3f[ot]*sqrt(FAo->result());
              break;
  }
  return(zi*zf);
}

Cplx matel::coeffS(void) const
{
  if(iff==VEC && off==VEC)
  {
    double mz = ival->get(MZ),
	   gz = ival->get(GamZ);
    return(vg0(it,*ival) * vg0(ot,*ival) / (mz*mz) * (1.+ I*gz/mz));
  }
  else
    return(Cplx(0));
}

Cplx matel::coeffSp(void) const
{
  if(iff==VEC && off==VEC)
  {
    double mz = ival->get(MZ);
    return(-vg0(it,*ival) * vg0(ot,*ival) / (mz*mz*mz*mz));
  }
  else
    return(Cplx(0));
}

Cplx matel::resoffZ(void) const
{
  double mz = ival->get(MZ);
  double gi = g0(it,iff,*ival), gf = g0(ot,off,*ival);
  return(gi*gf*sqr(1-s/(mz*mz))/s);
}

Cplx matel::result(void) const
{
  double mz = ival->get(MZ),
         gz = ival->get(GamZ);
  Cplx sminuss0(s - mz*mz, mz*gz);
  return(coeffR()/sminuss0 + coeffS() + coeffSp()*sminuss0 + resoffZ());
}

/*************************************************************************/

Cplx SW_SMLO::result(void) const
{
  return((1 - vz0(ftyp,*ival)/az0(ftyp,*ival))/(4*fabs(Qf[ftyp])));
}

Cplx FA_SMLO::result(void) const
{
  return(sqr(az0(ftyp,*ival)));
}

Cplx FV_SMLO::result(void) const
{
  double QVf = 1-4*fabs(Qf[ftyp])*realreg(sw->result());
  return(fa->result()*QVf*QVf);
}

Cplx FV_SMLO::errest(void) const
{
  double QVf0 = vz0(ftyp,*ival)/az0(ftyp,*ival),
         xsw0 = vz0(ftyp,*ival)*az0(ftyp,*ival)*fabs(Qf[ftyp]);
  return(sqrt(sqr(fa->errest()*QVf0*QVf0) + sqr(-8*xsw0 * sw->errest())));
}

/*************************************************************************/

extern const double Qf[5] = { -1., 0, +0.6666666666666667, 
                              -0.3333333333333333, -0.3333333333333333 };
extern const double I3f[5] = { -0.5, +0.5, +0.5, -0.5, -0.5 };

// tree-level axial-vector Z vertex factors
double az0(int type, const inval& input)
{
  double el = sqrt(4*Pi * input.get(al)),
         cw = input.get(MW)/input.get(MZ);
  double sw = sqrt(1-cw*cw);
  
  switch(type)
  {
    case LEP: 
    case DQU:
    case BQU: 
      return(-el/(4*sw*cw));
    case NEU:
    case UQU:
      return(+el/(4*sw*cw));
  }
  return 0;
}

// tree-level vector Z vertex factors
double vz0(int type, const inval& input)
{
  double cw = input.get(MW)/input.get(MZ);
  return az0(type,input) * (1- 4*fabs(Qf[type])*(1 - cw*cw));
}

double z0(int type, int formt, const inval& input)
{
  if(formt==VEC)
  { return vz0(type,input); }
  else
  { return az0(type,input); }
}

double g0(int type, int formt, const inval& input)
{
  if(formt==VEC)
  { return vg0(type,input); }
  else
  { return 0; }
}

/*************************************************************************/

double realreg(Cplx x)
{
  double r = x.real();
  if(isfinite(r))
    return(r);
  else
    return(0);
}

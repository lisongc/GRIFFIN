/*
  File: classes.cc
  Purpose: Implement typed-input tree-level SM form factors, pole coefficients
           and scattering amplitudes, retaining the original SM formulas.
  Merged-version authors: Jackson Wallace and Lisong Chen
  Last updated: 2026-11-16
*/

/*-----------------------------------------------------------------------------
classes.cc
Lisong Chen (lic114@pitt.edu), Ayres Freitas (afreitas@pitt.edu)
last revision: 15 Aug 2024
-------------------------------------------------------------------------------
basic classes for form factors and matrix elements, including SM LO predictions
-----------------------------------------------------------------------------*/

#include "griffin/core/classes.h"
#include "griffin/sm/ff0.h"

namespace griffin {

// the base classes compute results at LO:

Cplx matel::coeffR(void) const
{
  Cplx zi, zf;
  
  switch(iff)
  {
    case Current::vector: zi = 2*I3f[index(it)]*sqrt(FAi->result())
                       *(1.-4*fabs(Qf[index(it)])*realreg(SWi->result()));
              break;
    case Current::axial: zi = 2*I3f[index(it)]*sqrt(FAi->result());
              break;
    default:  zi = 0;
  }
  switch(off)
  {
    case Current::vector: zf = 2*I3f[index(ot)]*sqrt(FAo->result())
                       *(1.-4*fabs(Qf[index(ot)])*realreg(SWo->result()));
              break;
    case Current::axial: zf = 2*I3f[index(ot)]*sqrt(FAo->result());
              break;
    default:  zf = 0;
  }
  return(zi*zf);
}

Cplx matel::coeffS(void) const
{
  if(iff==Current::vector && off==Current::vector)
  {
    double mz = ival->get(InputPar::ZMassComplex),
	   gz = ival->get(InputPar::ZWidthComplex);
    return(vg0(it,*ival) * vg0(ot,*ival) / (mz*mz) * (1.+ I*gz/mz));
  }
  else
    return(Cplx(0));
}

Cplx matel::coeffSp(void) const
{
  if(iff==Current::vector && off==Current::vector)
  {
    double mz = ival->get(InputPar::ZMassComplex);
    return(-vg0(it,*ival) * vg0(ot,*ival) / (mz*mz*mz*mz));
  }
  else
    return(Cplx(0));
}

Cplx matel::resoffZ(void) const
{
  double mz = ival->get(InputPar::ZMassComplex);
  double gi = g0(it,iff,*ival), gf = g0(ot,off,*ival);
  double res = gi*gf/s;
  if(it==ot)
  {
    double t = -s/2*(1-cost);
    double g0v = g0(it,Current::vector,*ival), g0a = g0(it,Current::axial,*ival),
           z0v = z0(it,Current::vector,*ival), z0a = z0(it,Current::axial,*ival);
    if(iff==off)
    {
      if(iff==Current::vector || iff==Current::axial)  // iff=off=Current::vector or Current::axial
        res += ((z0v*z0v+z0a*z0a)/(t-mz*mz) + (g0v*g0v+g0a*g0a)/t)/2;
      else                      // iff=off=Current::scalar or Current::pseudoscalar
        res += ((z0v*z0v-z0a*z0a)/(t-mz*mz) + (g0v*g0v-g0a*g0a)/t);
    }
    else
    {
      if(iff==Current::vector || iff==Current::axial)  // iff=Current::vector and off=Current::axial, or vice versa
        res += (2*z0v*z0a/(t-mz*mz) + 2*g0v*g0a/t)/2;
    }
  }
  return(res);
}

Cplx matel::result(void) const
{
  double mz = ival->get(InputPar::ZMassComplex),
         gz = ival->get(InputPar::ZWidthComplex);
  Cplx sminuss0(s - mz*mz, mz*gz);
  return(coeffR()/sminuss0 + resoffZ());
}

/*************************************************************************/

Cplx SW_SMLO::result(void) const
{
  return((1 - vz0(ftyp,*ival)/az0(ftyp,*ival))/(4*fabs(Qf[index(ftyp)])));
}

Cplx FA_SMLO::result(void) const
{
  return(sqr(az0(ftyp,*ival)));
}

Cplx FV_SMLO::result(void) const
{
  double QVf = 1-4*fabs(Qf[index(ftyp)])*realreg(sw->result());
  return(fa->result()*QVf*QVf);
}

Cplx FV_SMLO::errest(void) const
{
  double QVf0 = vz0(ftyp,*ival)/az0(ftyp,*ival),
         xsw0 = vz0(ftyp,*ival)*az0(ftyp,*ival)*fabs(Qf[index(ftyp)]);
  return(sqrt(sqr(fa->errest()*QVf0*QVf0) + sqr(-8*xsw0 * sw->errest())));
}

/*************************************************************************/

// electric charges of different fermion types
extern const double Qf[20] = { 0, 
                               -0.3333333333333333,
			       +0.6666666666666667, 
                               -0.3333333333333333,
			       +0.6666666666666667, 
                               -0.3333333333333333,
			       +0.6666666666666667, 
                               -0.3333333333333333,
			       +0.6666666666666667, 
			       0,
			       0,
			       -1,
			       0,
			       -1,
			       0,
			       -1,
			       0,
			       -1,
			       0,
			       0 };

// weak isospin of different fermion types
extern const double I3f[20] = { 0,
				-0.5, 
				+0.5,
				-0.5, 
				+0.5,
				-0.5, 
				+0.5,
				-0.5, 
				+0.5,
				0,
				0,
				-0.5, 
				+0.5,
				-0.5, 
				+0.5,
				-0.5, 
				+0.5,
				-0.5, 
				+0.5,
				0 };
				

// tree-level axial-vector Z vertex factors
double az0(Fermion type, const inval& input)
{
  double el = sqrt(4*Pi * input.get(InputPar::alpha)),
         cw = input.get(InputPar::WMassComplex)/input.get(InputPar::ZMassComplex);
  double sw = sqrt(1-cw*cw);
  
  return(I3f[index(type)]*el/(2*sw*cw));
}

// tree-level vector Z vertex factors
double vz0(Fermion type, const inval& input)
{
  double cw = input.get(InputPar::WMassComplex)/input.get(InputPar::ZMassComplex);
  return az0(type,input) * (1- 4*fabs(Qf[index(type)])*(1 - cw*cw));
}

double z0(Fermion type, Current formt, const inval& input)
{
  if(formt==Current::vector)
  { return vz0(type,input); }
  if(formt==Current::axial)
  { return az0(type,input); }
  return 0;
}

double g0(Fermion type, Current formt, const inval& input)
{
  if(formt==Current::vector)
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

} // namespace

/*
  File: reference_scattering.h
  Purpose: Supply test-only scattering assemblies from existing SM kernels;
           this is not an installed API or an independent SM-kernel derivation.
  Authors: Jackson Wallace and Lisong Chen
  Last updated: 2026-11-16
*/

#ifndef GRIFFIN_TEST_REFERENCE_SCATTERING_H
#define GRIFFIN_TEST_REFERENCE_SCATTERING_H

// Test-only assembly from the existing low-level SM classes. This is not an
// alternate installed API or an independent validation of the SM kernels.
// The input masses come from the public finite solver and are separately
// checked against the unchanged Jackson archive by compare_aas.py.
#include "griffin/smeft/SMEFTCalculator.h"
#include "griffin/sm/EWPOZ2.h"
#include "griffin/sm/ff.h"
#include "griffin/sm/ff0.h"
#include "griffin/sm/xscaas.h"
#include "griffin/smeft/ff0_smeft.h"
#include "griffin/smeft/xsc_smeft.h"

#include <cmath>
#include <stdexcept>

namespace regression_reference {
using namespace griffin;
using namespace griffin::prediction;

enum class ScatteringContent { ArchivedSMNNLO, ExtendedSMaas };

struct PoleValues { Cplx fa, sw; };

inline PoleValues pole(Fermion f, const inval& input, const smeft::Input& coefficients)
{
  const double a = az0(index(f),input), v = vz0(index(f),input);
  const double da = smeft::directZAxialShift(f,input,coefficients);
  const double dv = smeft::directZVectorShift(f,input,coefficients);
  PoleValues result;
  result.fa = FA_SMNNLO(f,input).result()+2.0*a*da;
  result.sw = SW_SMNNLO(f,input).result()
      +(v*da/a-dv)/(4.0*std::fabs(Qf[index(f)])*a);
  return result;
}

inline Cplx observable(PseudoObservable which, Fermion f, const inval& input,
                       const smeft::Input& coefficients)
{
  const PoleValues p = pole(f,input,coefficients);
  if(which==PseudoObservable::FA) return p.fa;
  if(which==PseudoObservable::SW) return p.sw;
  if(which!=PseudoObservable::FV)
    throw std::invalid_argument("reference only implements FA/SW/FV");
  const int flavor=index(f);
  const double a=az0(flavor,input),v=vz0(flavor,input);
  const double imaginaryRatio=(a*(iz1f(flavor,VEC,input)+iz1b(flavor,VEC,input))
      -v*(iz1f(flavor,AXV,input)+iz1b(flavor,AXV,input)))/(a*a);
  const double realRatio=1.0-4.0*std::fabs(Qf[flavor])*real(p.sw);
  return p.fa*(realRatio*realRatio+imaginaryRatio*imaginaryRatio);
}

class FixedPole : public psobs {
  Cplx value_;
public:
  FixedPole(Cplx value,const inval& input) : psobs(input),value_(value) {}
  Cplx result() const override { return value_; }
};

inline SMEFTCalculator::Amplitudes contacts(const ScatteringPoint& point,
                                           const smeft::Input& coefficients)
{
  const Current ci[]={Current::vector,Current::axial,Current::vector,Current::axial};
  const Current cf[]={Current::vector,Current::vector,Current::axial,Current::axial};
  SMEFTCalculator::Amplitudes result;
  for(unsigned k=0;k<4;++k)
    result[k]=smeft::fourFermionContactAmplitude(point.initial,point.final,ci[k],cf[k],coefficients);
  return result;
}

inline SMEFTCalculator::Amplitudes amplitudes(const ScatteringPoint& point,
    const inval& input,const smeft::Input& coefficients,ScatteringContent content)
{
  const PoleValues initial=pole(point.initial,input,coefficients);
  const PoleValues final=pole(point.final,input,coefficients);
  const FixedPole fai(initial.fa,input),faf(final.fa,input);
  const FixedPole swi(initial.sw,input),swf(final.sw,input);
  const Current ci[]={Current::vector,Current::axial,Current::vector,Current::axial};
  const Current cf[]={Current::vector,Current::vector,Current::axial,Current::axial};
  SMEFTCalculator::Amplitudes result=contacts(point,coefficients);
  for(unsigned k=0;k<4;++k)
  {
    if(content==ScatteringContent::ArchivedSMNNLO)
      result[k]+=mat_SMNNLO(point.initial,point.final,ci[k],cf[k],
          fai,faf,swi,swf,point.s,point.cosTheta,input).result();
    else
      result[k]+=mat_SMaas(point.initial,point.final,ci[k],cf[k],
          fai,faf,swi,swf,point.s,point.cosTheta,input).result();
  }
  return result;
}

inline double crossSection(const ScatteringPoint& point,
                           const SMEFTCalculator::Amplitudes& amplitudes)
{
  double even=0.0;
  for(Cplx value:amplitudes) even+=std::norm(value);
  const double angular=(1.0+point.cosTheta*point.cosTheta)*even
      +4.0*point.cosTheta*real(amplitudes[0]*conj(amplitudes[3])
                            +amplitudes[2]*conj(amplitudes[1]));
  const double color=index(point.final)<6 ? 3.0 : 1.0;
  return color*point.s/(32.0*Pi)*point.unitConversion*angular;
}

} // namespace regression_reference
#endif

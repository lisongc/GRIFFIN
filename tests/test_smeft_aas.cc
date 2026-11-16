/*
  File: test_smeft_aas.cc
  Purpose: Check the finite calculator's mixed electroweak/QCD amplitude,
           full and SM results, and cross sections against a separate assembly.
  Authors: Jackson Wallace and Lisong Chen
  Last updated: 2026-11-16
*/

// Integration of the existing mixed electroweak/QCD SM amplitude into the
// finite SMEFT calculator. These checks independently assemble the public
// low-level ingredients; they are not an independent theory calculation.
#include "griffin/SMEFT.h"
#include "griffin/sm/ff0.h"
#include "griffin/sm/xscaas.h"
#include "griffin/smeft/ff0_smeft.h"
#include "griffin/smeft/xsc_smeft.h"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <string>

using namespace griffin;
using namespace griffin::prediction;

namespace {

typedef SMEFTCalculator::Amplitudes Amplitudes;
const smeft::Generation one=smeft::Generation::first;
const smeft::Generation two=smeft::Generation::second;
const smeft::Generation three=smeft::Generation::third;
const double pbPerGeV2=0.38937966e9;

void require(bool condition, const std::string& message)
{
  if(!condition) throw std::runtime_error(message);
}

void close(Cplx actual, Cplx expected, const std::string& message,
           double relative=3e-12, double absolute=1e-15)
{
  if(!std::isfinite(real(actual)) || !std::isfinite(imag(actual))
      || !std::isfinite(real(expected)) || !std::isfinite(imag(expected))
      || std::abs(actual-expected)
          > absolute+relative*std::max(std::abs(actual),std::abs(expected)))
    throw std::runtime_error(message);
}

SMval makeSM(double alphaS=0.1179)
{
  SMval input;
  input.set(InputPar::ZMassExperiment,91.1876);
  input.set(InputPar::WMassExperiment,80.377);
  input.set(InputPar::ZWidthExperiment,2.4952);
  input.set(InputPar::WWidthExperiment,2.085);
  input.set(InputPar::alpha,1.0/137.03599976);
  input.set(InputPar::alphaS,alphaS);
  input.set(InputPar::HiggsMass,125.1);
  input.set(InputPar::topMass,172.5);
  input.set(InputPar::bottomMass,2.87);
  input.set(InputPar::charmMass,0.652);
  input.set(InputPar::tauMass,1.777);
  input.set(InputPar::strangeMass,0.0);
  input.set(InputPar::upMass,0.0);
  input.set(InputPar::downMass,0.0);
  input.set(InputPar::muonMass,0.0);
  input.set(InputPar::electronMass,0.0);
  input.set(InputPar::DeltaAlpha,0.059);
  input.set(InputPar::fermiConstant,1.166379e-5);
  return input;
}

// 0: zero, 1: vertices, 2: contacts, 3: both. All supported final-state
// families receive nonzero direct coefficients in the appropriate cases.
smeft::Input coefficients(unsigned kind, double sign=1.0)
{
  smeft::Input result;
  result.setLambdaGeV(246.0);
  if(kind==1 || kind==3)
  {
    result.setCphiD(sign*.0003);
    result.setCphiWB(sign*.0002);
    for(smeft::Generation generation : {one,two,three})
    {
      result.setCphil1(generation,generation,sign*.0001);
      result.setCphil3(generation,generation,-sign*.00015);
      result.setCphie(generation,generation,sign*.00012);
      result.setCphiq1(generation,generation,sign*.00008);
      result.setCphiq3(generation,generation,-sign*.00011);
      result.setCphiu(generation,generation,sign*.00009);
      result.setCphid(generation,generation,-sign*.00007);
    }
  }
  if(kind==2 || kind==3)
  {
    for(smeft::Generation generation : {two,three})
    {
      result.setCll(one,one,generation,generation,sign*.0004);
      result.setCll(one,generation,generation,one,sign*.0001);
      result.setCle(one,one,generation,generation,-sign*.0002);
      result.setCle(generation,generation,one,one,sign*.0003);
      result.setCee(one,one,generation,generation,sign*.0002);
    }
    for(smeft::Generation generation : {one,two,three})
    {
      result.setClq1(one,one,generation,generation,sign*.0003);
      result.setClq3(one,one,generation,generation,-sign*.0001);
      result.setClu(one,one,generation,generation,sign*.0002);
      result.setCld(one,one,generation,generation,-sign*.0002);
      result.setCqe(generation,generation,one,one,sign*.00015);
      result.setCeu(one,one,generation,generation,-sign*.0001);
      result.setCed(one,one,generation,generation,sign*.0001);
    }
  }
  return result;
}

class FixedPoleInput : public psobs {
  Cplx value_;
public:
  FixedPoleInput(Cplx value, const inval& input) : psobs(input), value_(value) {}
  Cplx result() const override { return value_; }
};

Cplx explicitFA(Fermion f, const inval& input, const smeft::Input& eft)
{
  return FA_SMNNLO(f,input).result()
      +2.0*az0(index(f),input)*smeft::directZAxialShift(f,input,eft);
}

Cplx explicitSW(Fermion f, const inval& input, const smeft::Input& eft)
{
  const int type=index(f);
  const double a=az0(type,input), v=vz0(type,input);
  const double da=smeft::directZAxialShift(f,input,eft);
  const double dv=smeft::directZVectorShift(f,input,eft);
  return SW_SMNNLO(f,input).result()-(dv/a-v*da/(a*a))/(4.0*std::fabs(Qf[type]));
}

struct Assembly {
  Amplitudes total, previous, poleAddition, backgroundAddition, contacts;
};

Assembly assemble(const inval& input, const smeft::Input& eft,
                  const ScatteringPoint& point)
{
  const FixedPoleInput fai(explicitFA(point.initial,input,eft),input);
  const FixedPoleInput faf(explicitFA(point.final,input,eft),input);
  const FixedPoleInput swi(explicitSW(point.initial,input,eft),input);
  const FixedPoleInput swf(explicitSW(point.final,input,eft),input);
  const double mz=input.get(InputPar::ZMassComplex);
  const Cplx denominator(point.s-mz*mz,mz*input.get(InputPar::ZWidthComplex));
  const Current ci[]={Current::vector,Current::axial,Current::vector,Current::axial};
  const Current cf[]={Current::vector,Current::vector,Current::axial,Current::axial};
  Assembly result;
  for(unsigned k=0;k<4;++k)
  {
    // Use two concrete objects: a base reference to the extended object is
    // not a reliable way to request the old nonvirtual residue/remainder.
    const mat_SMNNLO oldMatrix(point.initial,point.final,ci[k],cf[k],
        fai,faf,swi,swf,point.s,point.cosTheta,input);
    const mat_SMaas extendedMatrix(point.initial,point.final,ci[k],cf[k],
        fai,faf,swi,swf,point.s,point.cosTheta,input);
    result.contacts[k]=smeft::fourFermionContactAmplitude(point.initial,point.final,ci[k],cf[k],eft);
    result.previous[k]=oldMatrix.result()+result.contacts[k];
    result.poleAddition[k]=(extendedMatrix.coeffR()-oldMatrix.coeffR())/denominator;
    result.backgroundAddition[k]=extendedMatrix.resoffZaas();
    result.total[k]=extendedMatrix.result()+result.contacts[k];
    close(result.total[k],result.previous[k]+result.poleAddition[k]
          +result.backgroundAddition[k],"explicit residue/background decomposition");
  }
  return result;
}

double square(const Amplitudes& a, const ScatteringPoint& point)
{
  double even=0.0;
  for(Cplx value:a) even+=std::norm(value);
  const bool quark=point.final==Fermion::u || point.final==Fermion::d
      || point.final==Fermion::s || point.final==Fermion::c || point.final==Fermion::b;
  const double c=point.cosTheta;
  return (quark?3.0:1.0)*point.s/(32.0*Pi)*point.unitConversion
      *((1.0+c*c)*even+4.0*c*real(a[0]*conj(a[3])+a[2]*conj(a[1])));
}

unsigned testScan(const ElectroweakSetup& setup)
{
  const smeft::Input zero=coefficients(0);
  double maxPoleAddition=0.0, maxBackgroundAddition=0.0, maxCrossSectionChange=0.0;
  unsigned points=0;
  for(unsigned kind=0;kind<=3;++kind)
    for(double sign : {-1.0,1.0})
    {
      if(kind==0 && sign<0.0) continue;
      const smeft::Input eft=coefficients(kind,sign);
      const SMEFTCalculator calculation(setup,eft);
      for(Fermion final : {Fermion::muon,Fermion::tau,Fermion::u,Fermion::d,
                            Fermion::s,Fermion::c,Fermion::b})
        for(double energy : {50.0,60.0,90.0,95.0,120.0})
          for(double cosine : {-.3,.5})
          {
            const ScatteringPoint point(Fermion::electron,final,energy*energy,cosine,pbPerGeV2);
            const Assembly full=assemble(calculation.resolvedInput(),eft,point);
            const Assembly sm=assemble(calculation.smInputValues(),zero,point);
            const auto actual=calculation.amplitudes(point);
            const auto actualSM=calculation.smAmplitudes(point);
            const auto actualContact=calculation.contactAmplitudes(point);
            for(unsigned k=0;k<4;++k)
            {
              close(actual[k],full.total[k],"public full amplitude does not use complete mat_SMaas");
              close(actualSM[k],sm.total[k],"public SM result does not use complete mat_SMaas");
              close(actualContact[k],full.contacts[k],"contact contribution changed in extended amplitude");
              maxPoleAddition=std::max(maxPoleAddition,std::abs(full.poleAddition[k]));
              maxBackgroundAddition=std::max(maxBackgroundAddition,std::abs(full.backgroundAddition[k]));
            }
            const auto sigma=calculation.differentialCrossSection(point);
            const double fullSquare=square(full.total,point), smSquare=square(sm.total,point);
            close(sigma.total(),fullSquare,"extended full angular square");
            close(sigma.smResult(),smSquare,"extended SM angular square");
            close(sigma.modelCorrection(),fullSquare-smSquare,
                  "extended finite total-minus-SM result",3e-11,2e-12);
            maxCrossSectionChange=std::max(maxCrossSectionChange,
                std::fabs(fullSquare-square(full.previous,point)));
            ++points;
          }
    }
  require(maxPoleAddition>1e-12,"scan did not exercise the extra imaginary pole residue");
  require(maxBackgroundAddition>1e-12,"scan did not exercise the extra off-resonance term");
  require(maxCrossSectionChange>1e-7,"scan has no measurable added SM contribution");
  return points;
}

void testZeroStrongCoupling()
{
  // Construct the SM input with alpha_s=0 before Gmu solving. Both the pole
  // phase and the running-coupling remainder must then reduce to mat_SMNNLO.
  const SMval direct=makeSM(0.0);
  const SMvalGmu derived(direct);
  const ElectroweakSetup alpha(EWInputScheme::AlphaMWMZ,direct);
  const ElectroweakSetup gmu(EWInputScheme::GmuAlphaMZ,derived);
  close(getalphas(50.0*50.0,91.1876*91.1876,0.0),0.0,
        "zero strong-coupling running limit",0.0,0.0);
  for(const ElectroweakSetup* setup : {&alpha,&gmu})
    for(unsigned kind : {0u,3u})
    {
      const smeft::Input eft=coefficients(kind);
      const SMEFTCalculator calculation(*setup,eft);
      for(Fermion final : {Fermion::muon,Fermion::b})
        for(double energy : {50.0,90.0,120.0})
        {
          const ScatteringPoint point(Fermion::electron,final,energy*energy,.5,pbPerGeV2);
          const Assembly expected=assemble(calculation.resolvedInput(),eft,point);
          const auto actual=calculation.amplitudes(point);
          for(unsigned k=0;k<4;++k)
          {
            close(expected.poleAddition[k],0.0,"alpha_s=0 left an added pole term",0.0,1e-22);
            close(expected.backgroundAddition[k],0.0,"alpha_s=0 left an added background term",0.0,0.0);
            close(actual[k],expected.previous[k],"alpha_s=0 failed to recover mat_SMNNLO");
          }
          close(calculation.differentialCrossSection(point).total(),square(expected.previous,point),
                "alpha_s=0 full cross section differs from old SM content");
        }
    }
}

} // namespace

int main()
{
  try
  {
    const SMval direct=makeSM();
    const SMvalGmu derived(direct);
    const ElectroweakSetup alpha(EWInputScheme::AlphaMWMZ,direct);
    const ElectroweakSetup gmu(EWInputScheme::GmuAlphaMZ,derived);
    const unsigned count=testScan(alpha)+testScan(gmu);
    testZeroStrongCoupling();
    std::cout << "Mixed electroweak/QCD SMEFT amplitude tests: PASS ("
              << count << " scan points plus zero-alpha_s limits)\n";
    return 0;
  }
  catch(const std::exception& error)
  {
    std::cerr << "Mixed electroweak/QCD SMEFT amplitude tests: FAIL: " << error.what() << '\n';
    return 1;
  }
}

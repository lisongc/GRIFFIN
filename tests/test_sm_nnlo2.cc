/*
  Test the v1.2 SM amplitude port against frozen, full-precision results from
  unmodified upstream ca2d523, and check its integration with the typed API.
  These are port/regression checks, not an independent theory validation.
*/
#include "griffin/sm/EWPOZ2.h"
#include "griffin/sm/SMval.h"
#include "griffin/sm/ff.h"
#include "griffin/sm/xscnnlo2.h"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <type_traits>

using namespace griffin;

namespace {

static_assert(std::is_constructible<mat_SMNNLO2,Fermion,Fermion,Current,Current,
    double,double,double,double,double,double,const inval&>::value,
    "mat_SMNNLO2 must inherit the typed numeric-form-factor constructor");
static_assert(!std::is_constructible<mat_SMNNLO2,int,int,int,int,
    double,double,double,double,double,double,const inval&>::value,
    "mat_SMNNLO2 must reject integer labels for numeric form factors");
static_assert(!std::is_constructible<mat_SMNNLO2,int,int,int,int,
    const psobs&,const psobs&,const psobs&,const psobs&,double,double,
    const inval&>::value,
    "mat_SMNNLO2 must reject integer labels for borrowed form factors");

struct ReferencePoint {
  int finalType, initialCurrent, finalCurrent;
  double energy, cosine;
  // Re/Im pairs: total, mat_SMaas total, two-fermion-loop remainder,
  // fermion-boson remainder, and photon-vacuum-polarization contribution.
  double values[10];
};

const ReferencePoint reference[] = {
#include "upstream_v12_reference/points.inc"
};

void require(bool condition, const std::string& message)
{
  if(!condition) throw std::runtime_error(message);
}

void close(Cplx actual, Cplx expected, const std::string& message,
           double relative=2e-10, double absolute=2e-17)
{
  if(!std::isfinite(real(actual)) || !std::isfinite(imag(actual))
      || !std::isfinite(real(expected)) || !std::isfinite(imag(expected))
      || std::abs(actual-expected)
          > absolute+relative*std::max(std::abs(actual),std::abs(expected)))
  {
    std::ostringstream details;
    details.precision(17);
    details << message << ": got " << actual << ", expected " << expected;
    throw std::runtime_error(details.str());
  }
}

SMval makeInput()
{
  SMval input;
  input.set(InputPar::ZMassExperiment,91.1876);
  input.set(InputPar::WMassExperiment,80.377);
  input.set(InputPar::ZWidthExperiment,2.4952);
  input.set(InputPar::WWidthExperiment,2.085);
  input.set(InputPar::alpha,1.0/137.03599976);
  input.set(InputPar::alphaS,0.1179);
  input.set(InputPar::HiggsMass,125.1);
  input.set(InputPar::topMass,172.5);
  input.set(InputPar::bottomMass,2.87);
  input.set(InputPar::charmMass,1.27);
  input.set(InputPar::tauMass,1.777);
  input.set(InputPar::strangeMass,0.342);
  input.set(InputPar::upMass,0.246);
  input.set(InputPar::downMass,0.246);
  input.set(InputPar::muonMass,0.1057);
  input.set(InputPar::electronMass,0.000511);
  input.set(InputPar::DeltaAlpha,0.059);
  return input;
}

void checkReference(const SMval& input)
{
  double largestAddition[3] = {0.0,0.0,0.0};
  for(const ReferencePoint& point : reference)
  {
    const Fermion final=static_cast<Fermion>(point.finalType);
    const Current ci=static_cast<Current>(point.initialCurrent);
    const Current cf=static_cast<Current>(point.finalCurrent);
    const FA_SMNNLO fai(Fermion::electron,input), faf(final,input);
    const SW_SMNNLO swi(Fermion::electron,input), swf(final,input);
    const double s=point.energy*point.energy;
    const mat_SMNNLO2 upgraded(Fermion::electron,final,ci,cf,
        fai,faf,swi,swf,s,point.cosine,input);
    const mat_SMaas previous(Fermion::electron,final,ci,cf,
        fai,faf,swi,swf,s,point.cosine,input);
    const Cplx values[] = {upgraded.result(),previous.result(),
        upgraded.resoffZ2f(),upgraded.resoffZfb(),upgraded.resvpol()};
    std::ostringstream label;
    label << "upstream v1.2 f=" << point.finalType << " currents="
          << point.initialCurrent << point.finalCurrent << " E="
          << point.energy << " cos=" << point.cosine;
    for(unsigned j=0;j<5;++j)
      close(values[j],Cplx(point.values[2*j],point.values[2*j+1]),
            label.str()+" component "+std::to_string(j));

    close(values[0],values[1]+values[2]+values[3]+values[4],
          label.str()+" total decomposition",3e-12,2e-18);
    close(upgraded.resoffZ(),previous.resoffZ()+values[2]+values[3],
          label.str()+" off-resonance decomposition",3e-12,2e-18);
    close(upgraded.coeffR(),previous.coeffR(),
          label.str()+" inherited pole residue",0.0,0.0);
    for(unsigned j=0;j<3;++j)
      largestAddition[j]=std::max(largestAddition[j],std::abs(values[j+2]));

    const double t=-s/2.0*(1.0-point.cosine);
    const Cplx defaults=B1s(Fermion::electron,final,ci,cf,s,t,input,1,1,1);
    const Cplx explicitFlags=B1s(Fermion::electron,final,ci,cf,s,t,input,1,1,1,1,1);
    close(defaults,explicitFlags,label.str()+" default box flags",0.0,0.0);
  }
  for(double largest : largestAddition)
    require(largest>1e-12,"reference points must exercise every new contribution");
}

void checkBorrowedObservableUpdate(const SMval& input)
{
  const FA_SMNNLO fa(Fermion::electron,input);
  const SW_SMNNLO sw(Fermion::electron,input);
  psobsfix initialFA(real(fa.result()),input);
  const psobsfix finalFA(real(fa.result()),input);
  const psobsfix initialSW(real(sw.result()),input), finalSW(real(sw.result()),input);
  const mat_SMNNLO2 borrowed(Fermion::electron,Fermion::muon,
      Current::vector,Current::vector,initialFA,finalFA,initialSW,finalSW,
      90.0*90.0,0.3,input);
  const Cplx before=borrowed.result();
  const double updated=real(fa.result())*1.01;
  initialFA.setval(updated);
  close(initialFA.result(),updated,"psobsfix setter",0.0,0.0);
  const psobsfix freshFA(updated,input);
  const mat_SMNNLO2 fresh(Fermion::electron,Fermion::muon,
      Current::vector,Current::vector,freshFA,finalFA,initialSW,finalSW,
      90.0*90.0,0.3,input);
  close(borrowed.result(),fresh.result(),"borrowed observable update",0.0,0.0);
  require(std::abs(before-borrowed.result())>1e-12,
          "updated borrowed form factor did not affect the amplitude");
}

} // namespace

int main()
{
  try
  {
    const SMval input=makeInput();
    checkReference(input);
    checkBorrowedObservableUpdate(input);
    std::cout << "SM NNLO2 upstream-reference and typed-API tests: PASS ("
              << sizeof(reference)/sizeof(reference[0]) << " points)\n";
    return 0;
  }
  catch(const std::exception& error)
  {
    std::cerr << "SM NNLO2 tests: FAIL: " << error.what() << '\n';
    return 1;
  }
}

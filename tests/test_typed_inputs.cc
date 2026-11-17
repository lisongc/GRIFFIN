/*
  File: test_typed_inputs.cc
  Purpose: Require typed-only input/particle/current APIs and verify that SM
           mass definitions and inline calculations do not depend on header order.
  Authors: Jackson Wallace and Lisong Chen
  Last updated: 2026-11-16
*/

// Compile this same test in both orders. Experimental and complex-pole labels
// must never change meaning when the mass-conversion headers are included.
#ifdef GRIFFIN_TEST_SMVAL_FIRST
#include "griffin/sm/SMvalG.h"
#include "griffin/sm/EWPOZ2.h"
#else
#include "griffin/sm/EWPOZ2.h"
#include "griffin/sm/SMvalG.h"
#endif
#include "griffin/SMEFT.h"
#include "griffin/sm/ff.h"
#include "griffin/sm/ff0.h"
#include "griffin/sm/tools.h"
#include "griffin/sm/xscaas.h"
#include "griffin/smeft/ff0_smeft.h"

#include <cmath>
#include <iostream>
#include <stdexcept>
#include <type_traits>
#include <utility>

#if defined(MW) || defined(MZ) || defined(MH) || defined(ME) || defined(MM) \
 || defined(ML) || defined(MD) || defined(MS) || defined(MB) || defined(MU) \
 || defined(MC) || defined(MT) || defined(al) || defined(als) \
 || defined(Delal) || defined(Delalhad) || defined(Gmu) || defined(GamW) \
 || defined(GamZ) || defined(MWc) || defined(MZc) || defined(GWc) || defined(GZc) \
 || defined(ELE) || defined(MUO) || defined(TAU) || defined(NUE) \
 || defined(NUM) || defined(NUT) || defined(DQU) || defined(UQU) \
 || defined(SQU) || defined(CQU) || defined(BQU) || defined(VEC) \
 || defined(AXV) || defined(SCA) || defined(PSC) \
 || defined(COMPPOLESCHEME) || defined(RUNWIDTHSCHEME)
#error Legacy input labels must not be defined by public headers
#endif

using namespace griffin;

namespace {

// Taking an unambiguous function address also rejects accidental reintroduction
// of an integer overload, even if the typed overload is still present.
static_assert(std::is_same<decltype(&inval::set),
    void (inval::*)(InputPar,double)>::value,"set must accept InputPar only");
static_assert(std::is_same<decltype(&inval::get),
    double (inval::*)(InputPar) const>::value,"get must accept InputPar only");
static_assert(std::is_same<decltype(&matel::setform),
    void (matel::*)(Current,Current)>::value,"setform must accept Current only");
static_assert(std::is_same<decltype(&matel::setftype),
    void (matel::*)(Fermion,Fermion)>::value,"matrix setftype must accept Fermion only");
static_assert(std::is_same<decltype(&SW_SMLO::setftype),
    void (SW_SMLO::*)(Fermion)>::value,"SW setftype must accept Fermion only");
static_assert(std::is_same<decltype(&FA_SMLO::setftype),
    void (FA_SMLO::*)(Fermion)>::value,"FA setftype must accept Fermion only");
static_assert(std::is_same<decltype(&FV_SMLO::setftype),
    void (FV_SMLO::*)(Fermion)>::value,"FV setftype must accept Fermion only");

typedef double (*Vertex)(Fermion,const inval&);
typedef double (*CurrentVertex)(Fermion,Current,const inval&);
typedef double (*DirectShift)(Fermion,const inval&,const smeft::Input&);
typedef Cplx (*Contact)(Fermion,Fermion,Current,Current,const smeft::Input&);
static_assert(std::is_same<decltype(&az0),Vertex>::value,"az0 must accept Fermion only");
static_assert(std::is_same<decltype(&vz0),Vertex>::value,"vz0 must accept Fermion only");
static_assert(std::is_same<decltype(&vg0),Vertex>::value,"vg0 must accept Fermion only");
static_assert(std::is_same<decltype(&z0),CurrentVertex>::value,"z0 must use typed labels");
static_assert(std::is_same<decltype(&g0),CurrentVertex>::value,"g0 must use typed labels");
static_assert(std::is_same<decltype(&iz1f),CurrentVertex>::value,"iz1f must use typed labels");
static_assert(std::is_same<decltype(&iz1b),CurrentVertex>::value,"iz1b must use typed labels");
static_assert(std::is_same<decltype(&smeft::directZAxialShift),DirectShift>::value,
    "direct axial shift must accept Fermion only");
static_assert(std::is_same<decltype(&smeft::directZVectorShift),DirectShift>::value,
    "direct vector shift must accept Fermion only");
static_assert(std::is_same<decltype(&smeft::fourFermionContactAmplitude),Contact>::value,
    "contacts must use typed labels only");
static_assert(std::is_same<decltype(&partzwidth),
    double (*)(FA_SMLO&,FV_SMLO&,Fermion,const inval&,WidthScheme)>::value,
    "partial width must use typed fermion and width-convention labels");
static_assert(std::is_same<decltype(&zwidth),
    double (*)(FA_SMLO&,FV_SMLO&,const inval&,WidthScheme)>::value,
    "total width must use a typed width-convention label");

template<class Observable> struct TypedObservable {
  static_assert(std::is_constructible<Observable,Fermion,const inval&>::value,
      "typed observable constructor missing");
  static_assert(!std::is_constructible<Observable,int,const inval&>::value,
      "integer observable constructor survived");
};
template<class Matrix> struct TypedMatrix {
  static_assert(std::is_constructible<Matrix,Fermion,Fermion,Current,Current,
      double,double,double,double,double,double,const inval&>::value,
      "typed fixed-input matrix constructor missing");
  static_assert(!std::is_constructible<Matrix,int,int,int,int,
      double,double,double,double,double,double,const inval&>::value,
      "integer fixed-input matrix constructor survived");
  static_assert(std::is_constructible<Matrix,Fermion,Fermion,Current,Current,
      const psobs&,const psobs&,const psobs&,const psobs&,double,double,const inval&>::value,
      "typed observable-input matrix constructor missing");
  static_assert(!std::is_constructible<Matrix,int,int,int,int,
      const psobs&,const psobs&,const psobs&,const psobs&,double,double,const inval&>::value,
      "integer observable-input matrix constructor survived");
};

// Instantiation, not merely typedefs, is required to evaluate the assertions.
TypedObservable<SW_SMLO> swLO;
TypedObservable<SW_SMNLO> swNLO;
TypedObservable<SW_SMNNLO> swNNLO;
TypedObservable<FA_SMLO> faLO;
TypedObservable<FA_SMNLO> faNLO;
TypedObservable<FA_SMNNLO> faNNLO;
TypedObservable<FV_SMLO> fvLO;
TypedObservable<FV_SMNLO> fvNLO;
TypedObservable<FV_SMNNLO> fvNNLO;
TypedMatrix<matel> loMatrix;
TypedMatrix<mat_SMNNLO> nnloMatrix;
TypedMatrix<mat_SMaas> aasMatrix;

static_assert(!std::is_convertible<int,InputPar>::value,"InputPar must be strongly typed");
static_assert(!std::is_convertible<int,Fermion>::value,"Fermion must be strongly typed");
static_assert(!std::is_convertible<int,Current>::value,"Current must be strongly typed");
static_assert(!std::is_convertible<int,WidthScheme>::value,"WidthScheme must be strongly typed");

void close(double actual,double expected,const char* message)
{
  if(!std::isfinite(actual) || !std::isfinite(expected)
      || std::fabs(actual-expected)>2e-15+2e-13*std::fabs(expected))
    throw std::runtime_error(message);
}

void massDefinitionsAndInlineResponse()
{
  SMval input;
  input.set(InputPar::ZMassExperiment,91.1876);
  input.set(InputPar::ZWidthExperiment,2.4952);
  input.set(InputPar::WMassExperiment,80.377);
  input.set(InputPar::WWidthExperiment,2.085);
  input.set(InputPar::alpha,1.0/137.03599976);
  input.set(InputPar::alphaS,0.1179);
  input.set(InputPar::HiggsMass,125.1);
  input.set(InputPar::topMass,172.5);
  input.set(InputPar::bottomMass,2.87);
  input.set(InputPar::DeltaAlpha,0.059);
  input.set(InputPar::fermiConstant,1.166379e-5);
  const double mw=80.377/std::sqrt(1.0+std::pow(2.085/80.377,2));
  const double mz=91.1876/std::sqrt(1.0+std::pow(2.4952/91.1876,2));
  close(input.get(InputPar::WMassExperiment),80.377,"experimental W mass slot changed");
  close(input.get(InputPar::ZMassExperiment),91.1876,"experimental Z mass slot changed");
  close(input.get(InputPar::WMassComplex),mw,"complex W mass conversion changed");
  close(input.get(InputPar::ZMassComplex),mz,"complex Z mass conversion changed");
  close(input.get(InputPar::WWidthComplex),2.085*mw/80.377,"complex W width changed");
  close(input.get(InputPar::ZWidthComplex),2.4952*mz/91.1876,"complex Z width changed");
  const double rho=delrho2a2(&input);
  close(SW_SMNNLO(Fermion::electron,input).drho2a2(),
      mw*mw/(mz*mz)*rho,"SW inline response depends on header order");
  close(dr_SMNNLO(input).drho2a2(),
      rho/(1.0-mz*mz/(mw*mw)),"Delta-r inline response depends on header order");
  FA_SMNNLO fa(Fermion::electron,input);
  SW_SMNNLO sw(Fermion::electron,input);
  fa.setftype(Fermion::muon);
  sw.setftype(Fermion::muon);
  close(real(fa.result()),real(FA_SMNNLO(Fermion::muon,input).result()),
      "typed FA fermion setter changed result");
  close(real(sw.result()),real(SW_SMNNLO(Fermion::muon,input).result()),
      "typed SW fermion setter changed result");
}
} // namespace

int main()
{
  try {
    massDefinitionsAndInlineResponse();
    std::cout << "Typed-only inputs and header-order independence: PASS\n";
    return 0;
  } catch(const std::exception& error) {
    std::cerr << "Typed-only input test failed: " << error.what() << '\n';
    return 1;
  }
}

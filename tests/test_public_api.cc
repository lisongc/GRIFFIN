/*
  File: test_public_api.cc
  Purpose: Test the public entry header, electroweak setup, local contacts and
           amplitude-to-cross-section conversion, including rejected inputs.
  Authors: Jackson Wallace and Lisong Chen
  Last updated: 2026-11-17
*/

// Public facade, input setup and full-square algebra tests. No generic engine.
#include "griffin/SMEFT.h"
#include "griffin/SMEFT.h" // the recommended facade must be safe to include twice
#include "griffin/sm/EWPOZ2.h"
#include "griffin/sm/deltar.h"
#include "griffin/sm/tools.h"
#include "griffin/sm/xscnnlo.h"
#include "griffin/smeft/xsc_smeft.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <utility>

using namespace griffin;
using namespace griffin::prediction;

namespace {

static_assert(std::is_same<decltype(std::declval<const FinitePrediction<double>&>().smResult()),
    const double&>::value,"SM result accessor must return the stored value");
static_assert(std::is_same<decltype(std::declval<const SMEFTCalculator&>().smInputValues()),
    const inval&>::value,"SM input accessor must return the stored input values");
static_assert(std::is_same<decltype(std::declval<const SMEFTCalculator&>().smAmplitudes(
    std::declval<const ScatteringPoint&>())),SMEFTCalculator::Amplitudes>::value,
    "SM amplitudes must preserve the public channel-array type");

void require(bool condition, const std::string& message)
{ if(!condition) throw std::runtime_error(message); }

void close(Cplx actual, Cplx expected, const std::string& message)
{
  if(!std::isfinite(real(actual)) || !std::isfinite(imag(actual))
      || std::abs(actual-expected)
          > 1.0e-15+3.0e-12*std::max(std::abs(actual),std::abs(expected)))
    throw std::runtime_error(message);
}

template<typename Callable>
void rejects(const Callable& call, const std::string& message)
{
  try { call(); }
  catch(const std::exception&) { return; }
  throw std::runtime_error(message);
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

void testSetupAndFacade()
{
  const SMval direct=makeInput();
  const SMvalGmu derived(direct);
  rejects([&]{ ElectroweakSetup bad(EWInputScheme::GmuAlphaMZ,direct); },
          "Gmu setup accepted directly supplied MW");
  rejects([&]{ ElectroweakSetup bad(EWInputScheme::AlphaMWMZ,derived); },
          "alpha setup accepted Gmu-derived inputs");
  const ElectroweakSetup alpha(EWInputScheme::AlphaMWMZ,direct);
  const ElectroweakSetup gmu(EWInputScheme::GmuAlphaMZ,derived);
  require(&alpha.input()==&direct && &gmu.input()==&derived,
          "setup does not describe its supplied input");
  require(alpha.inputScheme()==EWInputScheme::AlphaMWMZ
          && gmu.inputScheme()==EWInputScheme::GmuAlphaMZ,"setup scheme lost");
  smeft::Input zero;
  zero.setLambdaGeV(1000.0);
  const ScatteringPoint point(Fermion::electron,Fermion::muon,8100.0,0.5,0.38937966e9);
  require(point.initial==Fermion::electron && point.final==Fermion::muon,
          "scattering flavors lost");
  close(point.s,8100.0,"scattering invariant lost");
  close(point.cosTheta,0.5,"scattering cosine lost");
  close(point.unitConversion,0.38937966e9,"scattering conversion lost");
  const CurrentChannel vv(Current::vector,Current::vector);
  require(vv.initial==Current::vector && vv.final==Current::vector,
          "current-channel identity lost");
  for(const ElectroweakSetup* setup : {&alpha,&gmu})
  {
    // Inline SM responses must keep their pole-mass slots even when the
    // umbrella also makes the running-width input classes available.
    const inval& input=setup->input();
    const double mw=input.get(InputPar::WMassComplex);
    const double mz=input.get(InputPar::ZMassComplex);
    const double rho=delrho2a2(&input);
    close(SW_SMNNLO(Fermion::electron,input).drho2a2(),
          (mw*mw)/(mz*mz)*rho,"umbrella changed SW pole-mass prefactor");
    close(dr_SMNNLO(input).drho2a2(),
          rho/(1.0-(mz*mz)/(mw*mw)),"umbrella changed Delta-r pole-mass prefactor");
    const SMEFTCalculator calculation(*setup,zero);
    const auto sigma=calculation.differentialCrossSection(point);
    require(sigma.total()>0.0,"zero-coefficient cross section is not positive");
    close(sigma.total(),sigma.smResult(),"zero coefficients differ from SM result");
    close(sigma.modelCorrection(),0.0,"zero coefficients have a model change");
    close(calculation.amplitude(point,vv).total(),
          calculation.amplitudes(point)[0],"channel accessor differs from array");
    close(calculation.smInputValues().get(InputPar::WMassComplex),
          calculation.resolvedInput().get(InputPar::WMassComplex),
          "zero coefficients differ from SM input values");
    close(calculation.smAmplitudes(point)[0],calculation.amplitudes(point)[0],
          "zero coefficients differ from SM amplitudes");
  }
}

void testContactAmplitudeApi()
{
  smeft::Input coefficients;
  coefficients.setLambdaGeV(1000.0);
  const smeft::Generation one=smeft::Generation::first;
  const smeft::Generation two=smeft::Generation::second;
  coefficients.setCee(one,one,two,two,0.2);
  const double expected=0.25*0.2/(1000.0*1000.0);
  close(smeft::fourFermionContactAmplitude(Fermion::electron,Fermion::muon,
        Current::vector,Current::vector,coefficients),expected,
        "typed public contact amplitude normalization");
}

void testCrossSectionFromAmplitudes()
{
  typedef std::array<Cplx,4> Amplitudes;
  const double s=100.0,c=0.25,conversion=2.0;
  const auto leptonSquare=[&](const Amplitudes& values) {
    return differentialCrossSectionFromAmplitudes(
        values,Fermion::electron,Fermion::muon,s,c,conversion);
  };
  const Amplitudes a={{Cplx(1.0,0.2),Cplx(-0.4,0.3),Cplx(0.1,-0.5),Cplx(0.7,0.2)}};
  const Amplitudes b={{Cplx(0.2,-0.1),Cplx(0.3,0.1),Cplx(-0.2,0.4),Cplx(0.1,0.2)}};
  double even=0.0;
  for(Cplx value : a) even+=std::norm(value);
  const double expected=s/(32.0*Pi)*conversion
      *((1+c*c)*even+4*c*real(a[0]*conj(a[3])+a[2]*conj(a[1])));
  close(leptonSquare(a),expected,"angular-square normalization");
  for(Fermion final : {Fermion::muon,Fermion::tau})
    close(differentialCrossSectionFromAmplitudes(a,Fermion::electron,final,s,c,conversion),
          expected,"charged-lepton normalization");
  for(Fermion final : {Fermion::d,Fermion::u,Fermion::s,Fermion::c,Fermion::b})
    close(differentialCrossSectionFromAmplitudes(a,Fermion::electron,final,s,c,conversion),
          3.0*expected,"quark color multiplicity");
  close(differentialCrossSectionFromAmplitudes(a,Fermion::electron,Fermion::muon,s,c),
        expected/conversion,"default GeV^-2 normalization");
  close(leptonSquare(Amplitudes{}),0.0,"zero amplitudes gave nonzero cross section");
  for(std::size_t k=0; k<4; ++k)
  {
    Amplitudes single={};
    single[k]=Cplx(1.0,0.5);
    close(leptonSquare(single),s/(32.0*Pi)*conversion*(1+c*c)*1.25,
          "single VV/AV/VA/AA channel normalization");
  }
  // Compare the legacy angular expression with an independent sum of chiral
  // squares, then apply the two-body phase space and physical unit conversions.
  // Supplied finite amplitudes also support the angular endpoints.
  const Cplx vv=a[0],av=a[1],va=a[2],aa=a[3];
  const Cplx ll=vv+av+va+aa,rr=vv-av-va+aa;
  const Cplx lr=vv+av-va-aa,rl=vv-av+va-aa;
  for(double angle : {-1.0,-0.8,0.0,0.5,0.8,1.0})
  {
    const double angular=real((1+angle*angle)
        *(vv*conj(vv)+av*conj(av)+va*conj(va)+aa*conj(aa))
        +4*angle*(vv*conj(aa)+va*conj(av)));
    const double chiral=0.25*((1+angle)*(1+angle)*(std::norm(ll)+std::norm(rr))
        +(1-angle)*(1-angle)*(std::norm(lr)+std::norm(rl)));
    close(chiral,angular,"chiral-square identity");
    for(Fermion final : {Fermion::muon,Fermion::d})
    {
      const double colors=final==Fermion::d ? 3.0 : 1.0;
      const double sigma=differentialCrossSectionFromAmplitudes(
          a,Fermion::electron,final,s,angle);
      close(sigma*(32.0*Pi*s),colors*s*s*chiral,
            "spin-averaged, color-summed matrix-element normalization");
      for(double unit : {1.0,0.38937966e6,0.38937966e9})
        close(differentialCrossSectionFromAmplitudes(
              a,Fermion::electron,final,s,angle,unit),
              colors*s/(32.0*Pi)*angular*unit,"legacy GeV^-2/nb/pb normalization");
    }
  }
  Amplitudes plus=a,minus=a;
  for(std::size_t k=0; k<4; ++k) { plus[k]+=b[k]; minus[k]-=b[k]; }
  close(leptonSquare(plus)+leptonSquare(minus),
        2.0*(leptonSquare(a)+leptonSquare(b)),
        "full-square quadratic identity");
  require(leptonSquare(b)>0.0,"correction square was discarded");
  for(double bad : {0.0,-1.0})
  {
    rejects([&]{ differentialCrossSectionFromAmplitudes(a,Fermion::electron,
        Fermion::muon,bad,c); },"nonpositive s accepted");
    rejects([&]{ differentialCrossSectionFromAmplitudes(a,Fermion::electron,
        Fermion::muon,s,c,bad); },"nonpositive conversion accepted");
  }
  for(double bad : {-1.1,1.1})
    rejects([&]{ differentialCrossSectionFromAmplitudes(a,Fermion::electron,
        Fermion::muon,s,bad); },"invalid angle accepted");
  for(double bad : {std::numeric_limits<double>::quiet_NaN(),
                   std::numeric_limits<double>::infinity(),
                   -std::numeric_limits<double>::infinity()})
  {
    rejects([&]{ differentialCrossSectionFromAmplitudes(a,Fermion::electron,
        Fermion::muon,bad,c); },"nonfinite s accepted");
    rejects([&]{ differentialCrossSectionFromAmplitudes(a,Fermion::electron,
        Fermion::muon,s,bad); },"nonfinite angle accepted");
    rejects([&]{ differentialCrossSectionFromAmplitudes(a,Fermion::electron,
        Fermion::muon,s,c,bad); },"nonfinite unit conversion accepted");
  }
  for(Fermion bad : {Fermion::muon,Fermion::nuElectron,static_cast<Fermion>(999)})
    rejects([&]{ differentialCrossSectionFromAmplitudes(a,bad,Fermion::muon,s,c); },
            "unsupported initial state accepted");
  for(Fermion bad : {Fermion::electron,Fermion::nuElectron,Fermion::nuMuon,
                    Fermion::nuTau,static_cast<Fermion>(6),static_cast<Fermion>(999)})
    rejects([&]{ differentialCrossSectionFromAmplitudes(a,Fermion::electron,bad,s,c); },
            "unsupported final state accepted");
}

} // namespace

int main()
{
  try
  {
    testSetupAndFacade();
    testContactAmplitudeApi();
    testCrossSectionFromAmplitudes();
    std::cout << "Finite public API tests: PASS\n";
    return 0;
  }
  catch(const std::exception& error)
  {
    std::cerr << "Finite public API tests: FAIL: " << error.what() << '\n';
    return 1;
  }
}

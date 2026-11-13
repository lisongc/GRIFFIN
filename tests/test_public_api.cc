// Public facade, input setup and full-square algebra tests. No generic engine.
#include "griffin/SMEFT.h"
#include "griffin/SMEFT.h" // the recommended facade must be safe to include twice
#include "griffin/sm/EWPOZ2.h"
#include "griffin/sm/deltar.h"
#include "griffin/sm/xscnnlo.h"
#include "griffin/sm/xscnnlo.h"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <string>

using namespace griffin;
using namespace griffin::prediction;

namespace {

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
    close(sigma.total(),sigma.reference(),"zero coefficients differ from reference");
    close(sigma.modelCorrection(),0.0,"zero coefficients have a model change");
    close(calculation.amplitude(point,vv).total(),
          calculation.amplitudes(point)[0],"channel accessor differs from array");
  }
}

void testFermionPairCrossSection()
{
  typedef FermionPairCrossSection CrossSection;
  const double s=100.0,c=0.25,conversion=2.0;
  const CrossSection lepton(Fermion::electron,Fermion::muon,s,c,conversion);
  const CrossSection quark(Fermion::electron,Fermion::d,s,c,conversion);
  const Current ci[]={Current::vector,Current::axial,Current::vector,Current::axial};
  const Current cf[]={Current::vector,Current::vector,Current::axial,Current::axial};
  for(std::size_t k=0; k<4; ++k)
    require(CrossSection::initialCurrent(k)==ci[k] && CrossSection::finalCurrent(k)==cf[k],
            "current-channel ordering");
  CrossSection::Amplitudes a={{Cplx(1.0,0.2),Cplx(-0.4,0.3),Cplx(0.1,-0.5),Cplx(0.7,0.2)}};
  CrossSection::Amplitudes b={{Cplx(0.2,-0.1),Cplx(0.3,0.1),Cplx(-0.2,0.4),Cplx(0.1,0.2)}};
  double even=0.0;
  for(Cplx value : a) even+=std::norm(value);
  const double expected=s/(32.0*Pi)*conversion
      *((1+c*c)*even+4*c*real(a[0]*conj(a[3])+a[2]*conj(a[1])));
  close(lepton.squaredPrediction(a),expected,"angular-square normalization");
  close(quark.squaredPrediction(a),3.0*expected,"quark color multiplicity");
  CrossSection::Amplitudes plus=a,minus=a;
  for(std::size_t k=0; k<4; ++k) { plus[k]+=b[k]; minus[k]-=b[k]; }
  close(lepton.squaredPrediction(plus)+lepton.squaredPrediction(minus),
        2.0*(lepton.squaredPrediction(a)+lepton.squaredPrediction(b)),
        "full-square quadratic identity");
  require(lepton.squaredPrediction(b)>0.0,"correction square was discarded");
  rejects([]{ CrossSection::initialCurrent(4); },"invalid current-channel index accepted");
  rejects([]{ CrossSection bad(Fermion::electron,Fermion::muon,0.0,0.5); },"nonpositive s accepted");
  rejects([]{ CrossSection bad(Fermion::electron,Fermion::muon,100.0,1.1); },"invalid angle accepted");
  rejects([]{ CrossSection bad(Fermion::electron,Fermion::muon,100.0,0.5,0.0); },
          "nonpositive conversion accepted");
}

} // namespace

int main()
{
  try
  {
    testSetupAndFacade();
    testFermionPairCrossSection();
    std::cout << "Finite public API tests: PASS\n";
    return 0;
  }
  catch(const std::exception& error)
  {
    std::cerr << "Finite public API tests: FAIL: " << error.what() << '\n';
    return 1;
  }
}

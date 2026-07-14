#include <algorithm>
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <string>

#include "EWPOZ2.h"
#include "SMvalG.h"
#include "models/SMEFTProvider.h"
#include "models/SMProvider.h"
#include "prediction/PredictionEngine.h"
#include "smeft/EWPOZ_smeft.h"
#include "smeft/deltar_smeft.h"
#include "smeft/xsc_smeft.h"
#include "xscmassless.h"

using namespace griffin;
using namespace griffin::prediction;

namespace {

class ToyFullModelProvider : public TheoryProvider {
private:
  CalculationContext context_;

public:
  explicit ToyFullModelProvider(const CalculationContext& context)
    : context_(context) {}

  const CalculationContext& context() const override { return context_; }
  ProviderOutputKind outputKind() const override
  {
    return ProviderOutputKind::FullModel;
  }
  const char* name() const override { return "toy full model"; }

  PerturbativeSeries<Cplx> pseudoObservableSeries(
      PseudoObservable, Fermion) const override
  {
    PerturbativeSeries<Cplx> series;
    series.add(OrderTag::referenceSM(LoopOrder::LO),10.0,
               "toy full-model LO");
    series.add(OrderTag::referenceSM(LoopOrder::NLO),2.0,
               "toy full-model NLO increment");
    return series;
  }

  PerturbativeSeries<Cplx> amplitudeSeries(
      const ScatteringPoint&, const CurrentChannel&) const override
  {
    return PerturbativeSeries<Cplx>();
  }
};

void require(bool condition, const std::string& message)
{
  if(!condition)
    throw std::runtime_error(message);
}

void requireClose(Cplx actual, Cplx expected, const std::string& message,
                  double tolerance = 2.0e-12)
{
  const double scale = std::max(1.0,std::max(std::abs(actual),std::abs(expected)));
  if(std::abs(actual-expected) > tolerance*scale)
    throw std::runtime_error(message);
}

void requireClose(double actual, double expected, const std::string& message,
                  double tolerance = 2.0e-12)
{
  const double scale = std::max(1.0,std::max(std::fabs(actual),std::fabs(expected)));
  if(std::fabs(actual-expected) > tolerance*scale)
    throw std::runtime_error(message);
}

template<typename Callable>
std::string failureMessage(const Callable& callable,
                           const std::string& missingFailureMessage)
{
  try
  {
    callable();
  }
  catch(const std::exception& error)
  {
    return error.what();
  }
  throw std::runtime_error(missingFailureMessage);
}

SMval makeDirectInput()
{
  SMval input;
  input.set(MZ,91.1876);
  input.set(MW,80.377);
  input.set(GamZ,2.4952);
  input.set(GamW,2.085);
  input.set(al,1.0/137.03599976);
  input.set(als,0.1179);
  input.set(MH,125.1);
  input.set(MT,172.5);
  input.set(MB,2.87);
  input.set(MC,0.652);
  input.set(MS,0.0);
  input.set(MU,0.0);
  input.set(MD,0.0);
  input.set(ML,1.777);
  input.set(MM,0.0);
  input.set(ME,0.0);
  input.set(Delal,0.059);
  input.set(Gmu,1.166379e-5);
  return input;
}

smeft::Input makeEFT(double scale)
{
  const smeft::Generation one = smeft::Generation::first;
  const smeft::Generation two = smeft::Generation::second;
  smeft::Input input;
  input.setLambdaGeV(1000.0);
  input.setCphiD(0.02*scale);
  input.setCphiWB(-0.01*scale);
  input.setCphil1(one,one,0.03*scale);
  input.setCphil3(one,one,-0.02*scale);
  input.setCphie(one,one,0.01*scale);
  input.setCphil3(two,two,0.015*scale);
  input.setCll(one,two,two,one,0.004*scale);
  input.setCll(two,one,one,two,0.004*scale);
  input.setCll(one,one,two,two,0.012*scale);
  input.setCle(one,one,two,two,-0.008*scale);
  input.setCle(two,two,one,one,0.006*scale);
  input.setCee(one,one,two,two,0.005*scale);
  return input;
}

PredictionEngine makeEngine(const SMProvider& sm, const SMEFTProvider& eft)
{
  return PredictionEngine(sm,eft,CombinationRule::linearEFT());
}

void testMasslessCrossSectionBuilder()
{
  const double s = 100.0;
  const double cosine = 0.25;
  const double conversion = 2.0;
  const MasslessFermionCrossSectionBuilder leptonBuilder(
      Fermion::electron,Fermion::muon,s,cosine,conversion);
  const MasslessFermionCrossSectionBuilder quarkBuilder(
      Fermion::electron,Fermion::d,s,cosine,conversion);

  require(leptonBuilder.initialCurrent(0) == Current::vector
          && leptonBuilder.finalCurrent(0) == Current::vector
          && leptonBuilder.initialCurrent(1) == Current::axial
          && leptonBuilder.finalCurrent(1) == Current::vector
          && leptonBuilder.initialCurrent(2) == Current::vector
          && leptonBuilder.finalCurrent(2) == Current::axial
          && leptonBuilder.initialCurrent(3) == Current::axial
          && leptonBuilder.finalCurrent(3) == Current::axial,
          "massless current-channel ordering changed");

  MasslessFermionCrossSectionBuilder::Amplitudes reference = {};
  MasslessFermionCrossSectionBuilder::Amplitudes correction = {};
  reference[0] = 1.0;
  correction[0] = 2.0;
  const double prefactor = s/(32.0*Pi)*conversion;
  const double angularEven = 1.0+cosine*cosine;
  const double expectedSquare = prefactor*angularEven;
  const double expectedInterference = 4.0*expectedSquare;

  requireClose(leptonBuilder.squaredPrediction(reference),expectedSquare,
               "massless cross-section normalization changed");
  requireClose(leptonBuilder.linearInterference(reference,correction),
               expectedInterference,
               "massless linear interference normalization changed");
  requireClose(quarkBuilder.squaredPrediction(reference),3.0*expectedSquare,
               "massless quark color multiplicity changed");
}

void testPseudoObservables(const SMvalGmu& input,
                           const smeft::Input& eftInput)
{
  const CalculationContext context(EWInputScheme::GmuAlphaMZ,input);
  const SMProvider sm(context);
  const SMEFTProvider eft(context,eftInput);
  const PredictionEngine engine = makeEngine(sm,eft);
  const Fermion f = Fermion::electron;

  smeft::SW_SMEFTLO swLegacy(f,input,eftInput,EWInputScheme::GmuAlphaMZ);
  smeft::FA_SMEFTLO faLegacy(f,input,eftInput,EWInputScheme::GmuAlphaMZ);
  smeft::FV_SMEFTLO fvLegacy(f,input,eftInput,EWInputScheme::GmuAlphaMZ);

  const PredictionResult<Cplx> swLO = engine.pseudoObservable(
      PseudoObservable::SW,f,Accuracy::SMLO_plus_ModelLO());
  const PredictionResult<Cplx> faLO = engine.pseudoObservable(
      PseudoObservable::FA,f,Accuracy::SMLO_plus_ModelLO());
  const PredictionResult<Cplx> fvLO = engine.pseudoObservable(
      PseudoObservable::FV,f,Accuracy::SMLO_plus_ModelLO());

  requireClose(swLO.total(),swLegacy.result(),
               "new and legacy SMLO+SMEFTLO SW disagree");
  requireClose(faLO.total(),faLegacy.result(),
               "new and legacy SMLO+SMEFTLO FA disagree");
  requireClose(fvLO.total(),fvLegacy.result(),
               "new and legacy SMLO+SMEFTLO FV disagree");

  const PredictionResult<Cplx> swNNLO = engine.pseudoObservable(
      PseudoObservable::SW,f,Accuracy::SMNNLOPlus_plus_ModelLO());
  const PredictionResult<Cplx> faNNLO = engine.pseudoObservable(
      PseudoObservable::FA,f,Accuracy::SMNNLOPlus_plus_ModelLO());
  const PredictionResult<Cplx> fvNNLO = engine.pseudoObservable(
      PseudoObservable::FV,f,Accuracy::SMNNLOPlus_plus_ModelLO());

  requireClose(swNNLO.reference(),SW_SMNNLO(f,input).result(),
               "SM NNLOPlus SW adapter is inconsistent");
  requireClose(faNNLO.reference(),FA_SMNNLO(f,input).result(),
               "SM NNLOPlus FA adapter is inconsistent");
  requireClose(fvNNLO.reference(),FV_SMNNLO(f,input).result(),
               "SM NNLOPlus FV adapter is inconsistent");

  requireClose(swNNLO.modelCorrection(),swLO.modelCorrection(),
               "changing SM accuracy changed SMEFTLO SW");
  requireClose(faNNLO.modelCorrection(),faLO.modelCorrection(),
               "changing SM accuracy changed SMEFTLO FA");
  requireClose(fvNNLO.modelCorrection(),fvLO.modelCorrection(),
               "changing SM accuracy changed SMEFTLO FV");

  smeft::dr_SMEFTLO drLegacy(
      input,eftInput,EWInputScheme::GmuAlphaMZ);
  const PredictionResult<Cplx> dr = engine.pseudoObservable(
      PseudoObservable::DeltaR,f,Accuracy::SMNLO_plus_ModelLO());
  requireClose(dr.total(),drLegacy.result(),
               "new and legacy SMNLO+SMEFTLO Delta r disagree");

  const PredictionResult<Cplx> mw = engine.pseudoObservable(
      PseudoObservable::WMass,f,Accuracy::SMNNLOPlus_plus_ModelLO());
  requireClose(mw.reference(),input.get(InputPar::WMassComplex),
               "MW reference is not the context value");
  requireClose(mw.modelCorrection(),
               smeft::derivedMWShiftD6(input,eftInput).deltaMW,
               "MW D6 correction is inconsistent");
}

void testLinearity(const SMvalGmu& input)
{
  const CalculationContext context(EWInputScheme::GmuAlphaMZ,input);
  const SMProvider sm(context);
  const smeft::Input plusInput = makeEFT(+1.0);
  const smeft::Input minusInput = makeEFT(-1.0);
  const smeft::Input twiceInput = makeEFT(+2.0);
  const smeft::Input zeroInput = makeEFT(0.0);
  const SMEFTProvider plusEFT(context,plusInput);
  const SMEFTProvider minusEFT(context,minusInput);
  const SMEFTProvider twiceEFT(context,twiceInput);
  const SMEFTProvider zeroEFT(context,zeroInput);
  const PredictionEngine plus = makeEngine(sm,plusEFT);
  const PredictionEngine minus = makeEngine(sm,minusEFT);
  const PredictionEngine twice = makeEngine(sm,twiceEFT);
  const PredictionEngine zero = makeEngine(sm,zeroEFT);
  const Accuracy accuracy = Accuracy::SMNNLOPlus_plus_ModelLO();

  const PseudoObservable observables[] = {
    PseudoObservable::SW,
    PseudoObservable::FA,
    PseudoObservable::FV,
    PseudoObservable::DeltaR,
    PseudoObservable::WMass
  };
  for(unsigned int i = 0; i < sizeof(observables)/sizeof(observables[0]); ++i)
  {
    const Cplx p = plus.pseudoObservable(
        observables[i],Fermion::electron,accuracy).modelCorrection();
    const Cplx m = minus.pseudoObservable(
        observables[i],Fermion::electron,accuracy).modelCorrection();
    const Cplx t = twice.pseudoObservable(
        observables[i],Fermion::electron,accuracy).modelCorrection();
    const Cplx z = zero.pseudoObservable(
        observables[i],Fermion::electron,accuracy).modelCorrection();
    requireClose(z,0.0,
                 std::string("zero Wilson input failed for ")+
                 pseudoObservableName(observables[i]));
    requireClose(p+m,0.0,
                 std::string("sign reversal failed for ")+
                 pseudoObservableName(observables[i]));
    requireClose(t,2.0*p,
                 std::string("doubling failed for ")+
                 pseudoObservableName(observables[i]));
  }
}

void testCrossSection(const SMvalGmu& input,
                      const smeft::Input& eftInput)
{
  const CalculationContext context(EWInputScheme::GmuAlphaMZ,input);
  const SMProvider sm(context);
  const SMEFTProvider eft(context,eftInput);
  const PredictionEngine engine = makeEngine(sm,eft);
  const PredictionEngine smOnly(sm);
  const double gevToNb = 0.38937966e6;
  const ScatteringPoint point(Fermion::electron,Fermion::muon,
                              90.0*90.0,0.5,gevToNb);

  const PredictionResult<double> newLO = engine.differentialCrossSection(
      point,Accuracy::SMLO_plus_ModelLO());
  FA_SMLO faElectron(Fermion::electron,input);
  FA_SMLO faMuon(Fermion::muon,input);
  SW_SMLO swElectron(Fermion::electron,input);
  SW_SMLO swMuon(Fermion::muon,input);
  const smeft::DifferentialCrossSectionLinear legacy =
      smeft::differentialCrossSectionLinear(
          Fermion::electron,Fermion::muon,
          faElectron,faMuon,swElectron,swMuon,
          point.s,point.cosTheta,input,eftInput,
          EWInputScheme::GmuAlphaMZ,gevToNb);
  requireClose(newLO.reference(),legacy.sm,
               "new and legacy SMLO cross sections disagree");
  requireClose(newLO.modelCorrection(),legacy.d6,
               "new and legacy SMEFTLO cross-section corrections disagree");

  const PredictionResult<double> newNNLO = engine.differentialCrossSection(
      point,Accuracy::SMNNLOPlus_plus_ModelLO());
  requireClose(newNNLO.modelCorrection(),newLO.modelCorrection(),
               "changing SM accuracy changed strict SMEFTLO cross section");

  const smeft::Input minusInput = makeEFT(-1.0);
  const smeft::Input zeroInput = makeEFT(0.0);
  const SMEFTProvider minusEFT(context,minusInput);
  const SMEFTProvider zeroEFT(context,zeroInput);
  const PredictionEngine minus = makeEngine(sm,minusEFT);
  const PredictionEngine zero = makeEngine(sm,zeroEFT);
  const double minusD6 = minus.differentialCrossSection(
      point,Accuracy::SMNNLOPlus_plus_ModelLO()).modelCorrection();
  requireClose(newNNLO.modelCorrection()+minusD6,0.0,
               "linear cross section is not odd under C -> -C");
  requireClose(zero.differentialCrossSection(
                   point,Accuracy::SMNNLOPlus_plus_ModelLO()).modelCorrection(),
               0.0,"zero Wilson input changed the cross section");

  const ScatteringPoint unsupportedPoints[] = {
    ScatteringPoint(Fermion::electron,Fermion::electron,
                    90.0*90.0,0.5,gevToNb),
    ScatteringPoint(Fermion::electron,Fermion::nuElectron,
                    90.0*90.0,0.5,gevToNb),
    ScatteringPoint(Fermion::electron,Fermion::nuMuon,
                    90.0*90.0,0.5,gevToNb),
    ScatteringPoint(Fermion::muon,Fermion::tau,
                    90.0*90.0,0.5,gevToNb)
  };
  for(unsigned int i = 0;
      i < sizeof(unsupportedPoints)/sizeof(unsupportedPoints[0]); ++i)
  {
    const ScatteringPoint unsupported = unsupportedPoints[i];
    const std::string smFailure = failureMessage(
        [&]() {
          smOnly.differentialCrossSection(unsupported,Accuracy::SMLO());
        },"SM-only engine accepted an unsupported scattering process");
    const std::string modelFailure = failureMessage(
        [&]() {
          engine.differentialCrossSection(
              unsupported,Accuracy::SMLO_plus_ModelLO());
        },"SM+model engine accepted an unsupported scattering process");
    require(smFailure == modelFailure,
            "SM-only and SM+model process validation disagree");
  }

  const std::string legacyBhabhaFailure = failureMessage(
      [&]() {
        smeft::differentialCrossSectionLinear(
            Fermion::electron,Fermion::electron,
            faElectron,faElectron,swElectron,swElectron,
            point.s,point.cosTheta,input,eftInput,
            EWInputScheme::GmuAlphaMZ,gevToNb);
      },"legacy SMEFT cross section accepted Bhabha scattering");
  const std::string engineBhabhaFailure = failureMessage(
      [&]() {
        engine.differentialCrossSection(
            unsupportedPoints[0],Accuracy::SMLO_plus_ModelLO());
      },"prediction engine accepted Bhabha scattering");
  require(legacyBhabhaFailure == engineBhabhaFailure,
          "legacy and engine Bhabha validation disagree");
}

void testValidation(const SMval& directInput, const SMvalGmu& input,
                    const smeft::Input& eftInput)
{
  bool rejectedWrongGmuInput = false;
  try
  {
    CalculationContext invalid(EWInputScheme::GmuAlphaMZ,directInput);
  }
  catch(const std::invalid_argument&)
  {
    rejectedWrongGmuInput = true;
  }
  require(rejectedWrongGmuInput,
          "GmuAlphaMZ accepted an input that does not derive MW from Gmu");

  bool rejectedWrongDirectInput = false;
  try
  {
    CalculationContext invalid(EWInputScheme::AlphaMWMZ,input);
  }
  catch(const std::invalid_argument&)
  {
    rejectedWrongDirectInput = true;
  }
  require(rejectedWrongDirectInput,
          "AlphaMWMZ accepted an input that derives MW from Gmu");

  const CalculationContext gmu(EWInputScheme::GmuAlphaMZ,input);
  const CalculationContext direct(EWInputScheme::AlphaMWMZ,directInput);
  const SMProvider sm(gmu);
  const SMEFTProvider mismatched(direct,eftInput);

  bool rejectedContext = false;
  try
  {
    PredictionEngine invalid(sm,mismatched,CombinationRule::linearEFT());
  }
  catch(const std::invalid_argument&)
  {
    rejectedContext = true;
  }
  require(rejectedContext,"mismatched input schemes were not rejected");

  const SMEFTProvider eft(gmu,eftInput);
  const PredictionEngine engine = makeEngine(sm,eft);
  const PredictionEngine smOnly(sm);

  const std::string smSWFailure = failureMessage(
      [&]() {
        smOnly.pseudoObservable(
            PseudoObservable::SW,Fermion::nuElectron,Accuracy::SMLO());
      },"SM-only engine accepted neutrino SW");
  const std::string modelSWFailure = failureMessage(
      [&]() {
        engine.pseudoObservable(
            PseudoObservable::SW,Fermion::nuElectron,
            Accuracy::SMLO_plus_ModelLO());
      },"SM+model engine accepted neutrino SW");
  require(smSWFailure == modelSWFailure,
          "SM-only and SM+model neutrino SW validation disagree");
  const std::string legacySWFailure = failureMessage(
      [&]() {
        smeft::SW_SMEFTLO legacy(
            Fermion::nuElectron,input,eftInput,EWInputScheme::GmuAlphaMZ);
        legacy.result();
      },"legacy SMEFT observable accepted neutrino SW");
  require(smSWFailure == legacySWFailure,
          "legacy and prediction-engine neutrino SW validation disagree");

  bool rejectedOrder = false;
  try
  {
    engine.pseudoObservable(
        PseudoObservable::FA,Fermion::electron,
        Accuracy::smPlusModel(LoopOrder::NNLOPlus,LoopOrder::NLO,1));
  }
  catch(const std::domain_error&)
  {
    rejectedOrder = true;
  }
  require(rejectedOrder,"unavailable SMEFT NLO was not rejected");
}

void testFullModelDifference(const SMvalGmu& input)
{
  const CalculationContext context(EWInputScheme::GmuAlphaMZ,input);
  const SMProvider sm(context);
  const ToyFullModelProvider fullModel(context);
  const PredictionEngine engine(
      sm,fullModel,CombinationRule::fullModelDifference());
  const Accuracy accuracy = Accuracy::smPlusModel(
      LoopOrder::NNLOPlus,LoopOrder::NLO,0);
  const PredictionResult<Cplx> result = engine.pseudoObservable(
      PseudoObservable::FA,Fermion::electron,accuracy);
  const PerturbativeSeries<Cplx> smSeries = sm.pseudoObservableSeries(
      PseudoObservable::FA,Fermion::electron);
  const Cplx expected = smSeries.sumThrough(LoopOrder::NNLOPlus)
      + 12.0-smSeries.sumThrough(LoopOrder::NLO);
  requireClose(result.total(),expected,
               "full-model difference rule double counted the SM");
}

} // namespace

int main()
{
  try
  {
    const SMval directInput = makeDirectInput();
    const SMvalGmu input(directInput);
    const smeft::Input eftInput = makeEFT(+1.0);

    testMasslessCrossSectionBuilder();
    std::cout << "shared massless cross-section builder: PASS\n";
    testPseudoObservables(input,eftInput);
    std::cout << "legacy pseudo-observable and order separation: PASS\n";
    testLinearity(input);
    std::cout << "Wilson sign and doubling linearity: PASS\n";
    testCrossSection(input,eftInput);
    std::cout << "legacy/engine cross section and strict D6 truncation: PASS\n";
    testValidation(directInput,input,eftInput);
    std::cout << "input-scheme context and unavailable-order validation: PASS\n";
    testFullModelDifference(input);
    std::cout << "full-model difference combination rule: PASS\n";
    std::cout << "prediction architecture consistency checks: PASS\n";
    return 0;
  }
  catch(const std::exception& error)
  {
    std::cerr << "prediction architecture consistency checks: FAIL: "
              << error.what() << '\n';
    return 1;
  }
}

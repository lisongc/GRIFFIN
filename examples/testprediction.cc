/*
  Order-aware prediction-engine example.

  The same providers are reused for SM-only predictions and for three
  SM+SMEFT combinations. PredictionResult keeps the SM reference, the
  dimension-six correction, and their total separately.
*/

#include <iomanip>
#include <iostream>

#include "SMvalG.h"
#include "models/SMEFTProvider.h"
#include "models/SMProvider.h"
#include "prediction/PredictionEngine.h"

using namespace griffin;
using namespace griffin::prediction;

namespace {

SMval makeSMInput()
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
  input.set(Delal,0.059);
  input.set(Gmu,1.166379e-5);
  return input;
}

smeft::Input makeSMEFTInput()
{
  smeft::Input input;
  input.setLambdaGeV(1000.0);
  input.setCphiWB(0.01);
  input.setCphiD(-0.02);
  input.setCphil3(
      smeft::Generation::first,smeft::Generation::first,0.006);
  input.setCphil3(
      smeft::Generation::second,smeft::Generation::second,-0.010);
  return input;
}

void printFA(const PredictionEngine& prediction, const Accuracy& accuracy,
             const char* label)
{
  const PredictionResult<Cplx> result = prediction.pseudoObservable(
      PseudoObservable::FA,Fermion::muon,accuracy);
  std::cout << std::left << std::setw(30) << label
            << " SM=" << std::setw(20) << real(result.reference())
            << " D6=" << std::setw(20) << real(result.modelCorrection())
            << " total=" << real(result.total()) << '\n';
}

void printCrossSection(const PredictionEngine& prediction,
                       const ScatteringPoint& point,
                       const Accuracy& accuracy, const char* label)
{
  const PredictionResult<double> result =
      prediction.differentialCrossSection(point,accuracy);
  std::cout << std::left << std::setw(30) << label
            << " SM=" << std::setw(20) << result.reference()
            << " D6=" << std::setw(20) << result.modelCorrection()
            << " total=" << result.total() << '\n';
}

} // namespace

int main()
{
  const SMval directInput = makeSMInput();
  const SMvalGmu gmuInput(directInput);
  const smeft::Input eftInput = makeSMEFTInput();

  // One shared context enforces a common input scheme and SM reference point.
  const CalculationContext context(EWInputScheme::GmuAlphaMZ,gmuInput);
  const SMProvider sm(context);
  const SMEFTProvider eft(context,eftInput);
  const PredictionEngine smPrediction(sm);
  const PredictionEngine smeftPrediction(
      sm,eft,CombinationRule::linearEFT());

  std::cout << std::setprecision(15);
  std::cout << "FA(muon)\n";
  printFA(smPrediction,Accuracy::SMLO(),"SM LO");
  printFA(smPrediction,Accuracy::SMNLO(),"SM NLO");
  printFA(smPrediction,Accuracy::SMNNLOPlus(),"SM NNLOPlus");
  printFA(smeftPrediction,Accuracy::SMLO_plus_ModelLO(),
          "SMLO + SMEFTLO");
  printFA(smeftPrediction,Accuracy::SMNLO_plus_ModelLO(),
          "SMNLO + SMEFTLO");
  printFA(smeftPrediction,Accuracy::SMNNLOPlus_plus_ModelLO(),
          "SMNNLOPlus + SMEFTLO");

  const ScatteringPoint point(
      Fermion::electron,Fermion::muon,90.0*90.0,0.5,0.38937966e6);
  std::cout << "\nd sigma/d cos(theta) [nb]\n"
            << "e+e- -> mu+mu-, sqrt(s)=90 GeV, cos(theta)=0.5\n";
  printCrossSection(smPrediction,point,Accuracy::SMLO(),"SM LO");
  printCrossSection(smPrediction,point,Accuracy::SMNLO(),"SM NLO");
  printCrossSection(smPrediction,point,Accuracy::SMNNLOPlus(),
                    "SM NNLOPlus");
  printCrossSection(smeftPrediction,point,
                    Accuracy::SMLO_plus_ModelLO(),"SMLO + SMEFTLO");
  printCrossSection(smeftPrediction,point,
                    Accuracy::SMNLO_plus_ModelLO(),"SMNLO + SMEFTLO");
  printCrossSection(smeftPrediction,point,
                    Accuracy::SMNNLOPlus_plus_ModelLO(),
                    "SMNNLOPlus + SMEFTLO");
}

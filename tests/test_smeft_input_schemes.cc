/*
  Regression test for SMEFT input schemes.

  This file contains checks that are useful for developers but would obscure
  the public examples: constructor compatibility, linearization, the SM limit,
  Delta_mu isolation, an independent Delta-r formula, and the conversion
  identity between the two universal Z corrections.
*/

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>

#include "SMvalG.h"
#include "smeft/EWPOZ_smeft.h"
#include "smeft/deltar_smeft.h"
#include "smeft/ff0_smeft.h"

using namespace griffin;
using griffin::smeft::EWInputScheme;
using griffin::smeft::Generation;
using griffin::smeft::Input;

namespace {

void require(bool condition, const std::string& message)
{
  if(!condition)
    throw std::runtime_error(message);
}

bool close(Cplx left, Cplx right, double tolerance = 1e-14)
{
  return std::abs(left-right) <= tolerance;
}

SMval makeDirectMWInput()
{
  SMval input;
  input.set(MZ,91.1876);
  input.set(MW,80.358);
  input.set(al,1.0/137.03599976);
  input.set(als,0.1179);
  input.set(GamZ,2.4952);
  input.set(GamW,2.089);
  input.set(MH,125.1);
  input.set(MT,173.0);
  input.set(MB,2.87);
  input.set(Delal,0.05900);
  input.set(Gmu,1.166379e-5);
  return input;
}

Input makeSMEFTInput(double scale = 1.0)
{
  const Generation one = Generation::first;
  const Generation two = Generation::second;

  Input input;
  input.setLambdaGeV(1000.0);
  input.setCphiWB(scale*0.010);
  input.setCphiD(scale*-0.020);
  input.setCphil1(two,two,scale*0.015);
  input.setCphil3(two,two,scale*-0.010);
  input.setCphie(two,two,scale*0.008);
  input.setCphil3(one,one,scale*0.006);
  input.setCll(one,two,two,one,scale*0.004);
  input.setCll(two,one,one,two,scale*0.004);
  return input;
}

double expectedAlphaDeltaRD6(const inval& smInput, const Input& eftInput)
{
  const double alpha = smInput.get(InputPar::alpha);
  const double mw = smInput.get(InputPar::WMassComplex);
  const double mz = smInput.get(InputPar::ZMassComplex);
  const double cw = mw/mz;
  const double sw = std::sqrt(1.0-cw*cw);
  const double e = std::sqrt(4.0*Pi*alpha);
  const double vAlpha = 2.0*mw*sw/e;

  const Generation one = Generation::first;
  const Generation two = Generation::second;
  const double deltaMu = eftInput.Cphil3(one,one)
      + eftInput.Cphil3(two,two)
      - 0.5*(eftInput.Cll(one,two,two,one)
             + eftInput.Cll(two,one,one,two));
  const double bracket = 2.0*cw/sw*eftInput.CphiWB()
      + cw*cw/(2.0*sw*sw)*eftInput.CphiD() + deltaMu;
  return vAlpha*vAlpha*eftInput.inverseLambdaSquaredGeV()*bracket;
}

double smMWFixedPointMap(const inval& smInput, double mwSquared)
{
  inval shifted(smInput);
  shifted.set(InputPar::WMassComplex,std::sqrt(mwSquared));
  const double alpha = smInput.get(InputPar::alpha);
  const double mz = smInput.get(InputPar::ZMassComplex);
  const double gmu = smInput.get(InputPar::fermiConstant);
  const double k = Pi*alpha/(std::sqrt(2.0)*gmu*mz*mz);
  const double deltaRSM = real(dr_SMNNLO(shifted).result());
  return mz*mz*(0.5+std::sqrt(0.25-k*(1.0+deltaRSM)));
}

double solveShiftedMWSquared(const inval& smInput, double inputDeltaX)
{
  const double mz = smInput.get(InputPar::ZMassComplex);
  double value = std::pow(smInput.get(InputPar::WMassComplex),2);
  for(int iteration = 0; iteration < 100; ++iteration)
  {
    const double next = smMWFixedPointMap(smInput,value)
        + mz*mz*inputDeltaX;
    if(std::fabs(next-value) < 1.0e-9)
      return next;
    value = next;
  }
  throw std::runtime_error("shifted-MW fixed point did not converge");
}

} // namespace

int main()
{
  try
  {
    const SMval alphaMWInput = makeDirectMWInput();
    const SMvalGmu gmuInput(alphaMWInput);
    const Input eftInput = makeSMEFTInput();
    const Input oppositeEFTInput = makeSMEFTInput(-1.0);
    const Fermion fermion = Fermion::muon;
    const Generation one = Generation::first;
    const Generation two = Generation::second;

    SW_SMLO swSMGmu(fermion,gmuInput);
    SW_SMLO swSMAlpha(fermion,alphaMWInput);
    FA_SMLO faSMGmu(fermion,gmuInput);
    FV_SMLO fvSMGmu(fermion,gmuInput);

    smeft::SW_SMEFTLO swEFTGmu(
        fermion,gmuInput,eftInput,EWInputScheme::GmuAlphaMZ);
    smeft::SW_SMEFTLO swEFTAlpha(
        fermion,alphaMWInput,eftInput,EWInputScheme::AlphaMWMZ);
    smeft::FA_SMEFTLO faEFTGmu(
        fermion,gmuInput,eftInput,EWInputScheme::GmuAlphaMZ);
    smeft::FV_SMEFTLO fvEFTGmu(
        fermion,gmuInput,eftInput,EWInputScheme::GmuAlphaMZ);

    require(close(swEFTGmu.result()-swSMGmu.result(),
                  swEFTGmu.resD6Tree()),
            "GmuAlphaMZ SW failed linearization");
    require(close(swEFTAlpha.result()-swSMAlpha.result(),
                  swEFTAlpha.resD6Tree()),
            "AlphaMWMZ SW failed linearization");
    require(close(faEFTGmu.result()-faSMGmu.result(),
                  faEFTGmu.resD6Tree()),
            "FA failed linearization");
    require(close(fvEFTGmu.result()-fvSMGmu.result(),
                  fvEFTGmu.resD6Tree()),
            "FV failed linearization");

    // Exercise both signed currents for every supported fermion family and
    // generation.  Opposite Wilson coefficients must produce the exact
    // opposite result; no Wilson-dependent MW is inserted nonlinearly.
    const Fermion allFermions[] = {
        Fermion::d,Fermion::u,Fermion::s,Fermion::c,Fermion::b,
        Fermion::electron,Fermion::nuElectron,Fermion::muon,
        Fermion::nuMuon,Fermion::tau,Fermion::nuTau};
    for(const Fermion currentFermion : allFermions)
    {
      const double axial = smeft::azD6Tree(
          currentFermion,gmuInput,eftInput,EWInputScheme::GmuAlphaMZ);
      const double axialOpposite = smeft::azD6Tree(
          currentFermion,gmuInput,oppositeEFTInput,
          EWInputScheme::GmuAlphaMZ);
      const double vector = smeft::vzD6Tree(
          currentFermion,gmuInput,eftInput,EWInputScheme::GmuAlphaMZ);
      const double vectorOpposite = smeft::vzD6Tree(
          currentFermion,gmuInput,oppositeEFTInput,
          EWInputScheme::GmuAlphaMZ);
      require(std::isfinite(axial) && std::isfinite(vector),
              "non-finite SMEFT current");
      require(close(axial,-axialOpposite)
                  && close(vector,-vectorOpposite),
              "SMEFT current is not exactly linear");
    }

    // Check the two retained constructor APIs and the default scheme.
    smeft::SW_SMEFTLO swMacroDefault(MUO,gmuInput,eftInput);
    smeft::SW_SMEFTLO swEnumDefault(fermion,gmuInput,eftInput);
    smeft::SW_SMEFTLO swMacroAlpha(
        MUO,alphaMWInput,eftInput,EWInputScheme::AlphaMWMZ);
    require(close(swMacroDefault.result(),swEFTGmu.result()),
            "macro constructor default is not GmuAlphaMZ");
    require(close(swEnumDefault.result(),swEFTGmu.result()),
            "enum constructor default is not GmuAlphaMZ");
    require(close(swMacroAlpha.result(),swEFTAlpha.result()),
            "macro and enum AlphaMWMZ constructors disagree");

    // Delta_mu affects the GmuAlphaMZ current through the derived-MW shift,
    // but it is absent when MW is an independent AlphaMWMZ input.
    Input noCllInput = eftInput;
    noCllInput.setCll(one,two,two,one,0.0);
    noCllInput.setCll(two,one,one,two,0.0);
    smeft::SW_SMEFTLO swAlphaNoCll(
        fermion,alphaMWInput,noCllInput,EWInputScheme::AlphaMWMZ);
    smeft::SW_SMEFTLO swGmuNoCll(
        fermion,gmuInput,noCllInput,EWInputScheme::GmuAlphaMZ);
    require(close(swAlphaNoCll.resD6Tree(),swEFTAlpha.resD6Tree()),
            "Delta_mu leaked into AlphaMWMZ Z factors");
    require(!close(swGmuNoCll.resD6Tree(),swEFTGmu.resD6Tree()),
            "Delta_mu is absent from the GmuAlphaMZ derived-MW shift");

    // Zero Wilson coefficients recover the matching SM prediction.
    Input zeroInput;
    zeroInput.setLambdaGeV(1000.0);
    smeft::SW_SMEFTLO swZeroGmu(
        fermion,gmuInput,zeroInput,EWInputScheme::GmuAlphaMZ);
    smeft::SW_SMEFTLO swZeroAlpha(
        fermion,alphaMWInput,zeroInput,EWInputScheme::AlphaMWMZ);
    smeft::FA_SMEFTLO faZeroGmu(
        fermion,gmuInput,zeroInput,EWInputScheme::GmuAlphaMZ);
    smeft::FV_SMEFTLO fvZeroGmu(
        fermion,gmuInput,zeroInput,EWInputScheme::GmuAlphaMZ);
    require(close(swZeroGmu.result(),swSMGmu.result()),
            "zero coefficients failed GmuAlphaMZ SW limit");
    require(close(swZeroAlpha.result(),swSMAlpha.result()),
            "zero coefficients failed AlphaMWMZ SW limit");
    require(close(faZeroGmu.result(),faSMGmu.result()),
            "zero coefficients failed FA limit");
    require(close(fvZeroGmu.result(),fvSMGmu.result()),
            "zero coefficients failed FV limit");

    dr_SMNLO drSM(alphaMWInput);
    smeft::dr_SMEFTLO drEFT(
        alphaMWInput,eftInput,EWInputScheme::AlphaMWMZ);
    require(close(drEFT.result()-drSM.result(),drEFT.resD6Tree()),
            "Delta r failed linearization");
    require(close(drEFT.resD6Tree(),
                  Cplx(expectedAlphaDeltaRD6(alphaMWInput,eftInput))),
            "Delta-r formula failed its independent analytic check");

    // Compare the explicit linear response with an independent central
    // difference of the full Jackson-style fixed-point equation.  The latter
    // is used only as a regression oracle; observables use the analytic
    // linear displacement and never iterate with Wilson coefficients.
    const smeft::DerivedMWShiftD6 mwShift =
        smeft::derivedMWShiftD6(gmuInput,eftInput);
    const double plusMWSquared = solveShiftedMWSquared(
        gmuInput,+mwShift.inputDeltaMW2OverMZ2);
    const double minusMWSquared = solveShiftedMWSquared(
        gmuInput,-mwShift.inputDeltaMW2OverMZ2);
    const double fixedPointLinearShift =
        0.5*(plusMWSquared-minusMWSquared);
    require(std::fabs(fixedPointLinearShift-mwShift.deltaMW2)
                <= 1.0e-6*std::max(1.0,std::fabs(mwShift.deltaMW2)),
            "linear derived-MW response disagrees with fixed-point solve");

    std::cout << std::setprecision(16);
    std::cout << "Internal SMEFT input-scheme regression test\n";
    std::cout << "===========================================\n";
    std::cout << "Linearization checks        = PASS\n";
    std::cout << "All-fermion current checks  = PASS\n";
    std::cout << "Constructor checks          = PASS\n";
    std::cout << "SM-limit checks             = PASS\n";
    std::cout << "Delta-mu isolation checks   = PASS\n";
    std::cout << "Delta-r analytic check      = PASS\n";
    std::cout << "Raw delta(MW^2/MZ^2)        = "
              << mwShift.inputDeltaMW2OverMZ2 << '\n';
    std::cout << "SM fixed-point response     = "
              << mwShift.responseFactor << '\n';
    std::cout << "Linear delta(MW^2)          = "
              << mwShift.deltaMW2 << '\n';
    std::cout << "Fixed-point central shift   = "
              << fixedPointLinearShift << '\n';
    std::cout << "Derived-MW response check   = PASS\n";
    std::cout << "Overall result              = PASS\n";
  }
  catch(const std::exception& error)
  {
    std::cerr << "__testSMEFTSchemes failed: " << error.what() << '\n';
    return 1;
  }

  return 0;
}

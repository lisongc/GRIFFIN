/*
  User-facing example: SMEFT Z-pole form factors in two EW input schemes.

  This example intentionally contains no internal regression machinery.  It
  shows how to build the two SM inputs, assign Wilson coefficients, construct
  the SMEFT pseudo-observables, and compare their numerical predictions.
*/

#include <cmath>
#include <iomanip>
#include <iostream>

#include "SMvalG.h"
#include "smeft/EWPOZ_smeft.h"
#include "smeft/ff0_smeft.h"

using namespace griffin;
using griffin::smeft::EWInputScheme;
using griffin::smeft::Generation;
using griffin::smeft::Input;

namespace {

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

Input makeSMEFTInput()
{
  const Generation one = Generation::first;
  const Generation two = Generation::second;

  Input input;
  input.setLambdaGeV(1000.0);
  input.setCphiWB(0.010);
  input.setCphiD(-0.020);
  input.setCphil1(two,two,0.015);
  input.setCphil3(two,two,-0.010);
  input.setCphie(two,two,0.008);
  input.setCphil3(one,one,0.006);
  input.setCll(one,two,two,one,0.004);
  input.setCll(two,one,one,two,0.004);
  return input;
}

} // namespace

int main()
{
  // AlphaMWMZ uses the directly supplied MW.  GmuAlphaMZ copies the common
  // inputs and replaces MW by the SM value derived from Gmu.
  const SMval alphaMWInput = makeDirectMWInput();
  const SMvalGmu gmuInput(alphaMWInput);
  const Input eftInput = makeSMEFTInput();
  const smeft::DerivedMWShiftD6 derivedMWShift =
      smeft::derivedMWShiftD6(gmuInput,eftInput);
  const Fermion fermion = Fermion::muon;

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

  const double absoluteSWDifference =
      std::abs(swEFTGmu.result()-swEFTAlpha.result());
  const double absoluteSWD6Difference =
      std::abs(swEFTGmu.resD6Tree()-swEFTAlpha.resD6Tree());

  std::cout << std::setprecision(16);
  std::cout << "SM + linear SMEFT Z-pole example\n";
  std::cout << "================================\n";
  std::cout << "Fermion                       = muon (generation 2)\n";
  std::cout << "Lambda [GeV]                  = "
            << eftInput.lambdaGeV() << "\n\n";

  std::cout << "GmuAlphaMZ scheme             = "
            << smeft::inputSchemeName(EWInputScheme::GmuAlphaMZ) << '\n';
  std::cout << "Derived MW pole [GeV]         = "
            << gmuInput.get(InputPar::WMassComplex) << '\n';
  std::cout << "Linear SMEFT delta MW [GeV]   = "
            << derivedMWShift.deltaMW << '\n';
  std::cout << "SM + linear SMEFT MW [GeV]    = "
            << derivedMWShift.linearizedMW() << '\n';
  std::cout << "SM fixed-point response       = "
            << derivedMWShift.responseFactor << '\n';
  std::cout << "SW_SMLO                       = " << swSMGmu.result() << '\n';
  std::cout << "SW D6 tree                    = "
            << swEFTGmu.resD6Tree() << '\n';
  std::cout << "SW_SMLO + SMEFTLO             = "
            << swEFTGmu.result() << "\n\n";

  std::cout << "AlphaMWMZ scheme              = "
            << smeft::inputSchemeName(EWInputScheme::AlphaMWMZ) << '\n';
  std::cout << "Direct MW pole [GeV]          = "
            << alphaMWInput.get(InputPar::WMassComplex) << '\n';
  std::cout << "SW_SMLO                       = " << swSMAlpha.result() << '\n';
  std::cout << "SW D6 tree                    = "
            << swEFTAlpha.resD6Tree() << '\n';
  std::cout << "SW_SMLO + SMEFTLO             = "
            << swEFTAlpha.result() << "\n\n";

  std::cout << "Absolute full-SW scheme difference = "
            << absoluteSWDifference << '\n';
  std::cout << "Absolute D6-SW scheme difference   = "
            << absoluteSWD6Difference << "\n\n";

  std::cout << "GmuAlphaMZ azD6Tree           = "
            << smeft::azD6Tree(fermion,gmuInput,eftInput,
                                EWInputScheme::GmuAlphaMZ) << '\n';
  std::cout << "GmuAlphaMZ vzD6Tree           = "
            << smeft::vzD6Tree(fermion,gmuInput,eftInput,
                                EWInputScheme::GmuAlphaMZ) << '\n';
  std::cout << "AlphaMWMZ azD6Tree            = "
            << smeft::azD6Tree(fermion,alphaMWInput,eftInput,
                                EWInputScheme::AlphaMWMZ) << '\n';
  std::cout << "AlphaMWMZ vzD6Tree            = "
            << smeft::vzD6Tree(fermion,alphaMWInput,eftInput,
                                EWInputScheme::AlphaMWMZ) << "\n\n";

  std::cout << "FA_SMLO (GmuAlphaMZ)          = " << faSMGmu.result() << '\n';
  std::cout << "FA D6 tree                    = "
            << faEFTGmu.resD6Tree() << '\n';
  std::cout << "FA_SMLO + SMEFTLO             = "
            << faEFTGmu.result() << "\n\n";

  std::cout << "FV_SMLO (GmuAlphaMZ)          = " << fvSMGmu.result() << '\n';
  std::cout << "FV D6 tree                    = "
            << fvEFTGmu.resD6Tree() << '\n';
  std::cout << "FV_SMLO + SMEFTLO             = "
            << fvEFTGmu.result() << '\n';

  return 0;
}

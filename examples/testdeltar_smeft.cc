/*
  User-facing example: Delta r in the direct {alpha, MW, MZ} input scheme.

  In this scheme MW is an input, while Delta r (or equivalently the associated
  Gmu relation) is a prediction.  dr_SMEFTLO::result() returns SM NLO plus the
  linear dimension-six tree contribution; resD6Tree() returns that correction
  by itself.
*/

#include <cmath>
#include <iomanip>
#include <iostream>

#include "SMval.h"
#include "smeft/deltar_smeft.h"

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
  input.setCphil3(one,one,0.006);
  input.setCphil3(two,two,-0.010);
  input.setCll(one,two,two,one,0.004);
  input.setCll(two,one,one,two,0.004);
  return input;
}

} // namespace

int main()
{
  const SMval smInput = makeDirectMWInput();
  const Input eftInput = makeSMEFTInput();

  dr_SMNLO drSM(smInput);
  smeft::dr_SMEFTLO drEFT(
      smInput,eftInput,EWInputScheme::AlphaMWMZ);

  const double absoluteDifference =
      std::abs(drEFT.result()-drSM.result());
  const double relativeDifference =
      absoluteDifference/std::abs(drSM.result());

  std::cout << std::setprecision(16);
  std::cout << "SM + linear SMEFT Delta-r example\n";
  std::cout << "=================================\n";
  std::cout << "Input scheme                  = "
            << smeft::inputSchemeName(EWInputScheme::AlphaMWMZ) << '\n';
  std::cout << "Lambda [GeV]                  = "
            << eftInput.lambdaGeV() << '\n';
  std::cout << "Direct MW pole [GeV]          = "
            << smInput.get(InputPar::WMassComplex) << "\n\n";

  std::cout << "dr_SMNLO                      = " << drSM.result() << '\n';
  std::cout << "dr D6 tree                    = " << drEFT.resD6Tree() << '\n';
  std::cout << "dr_SMNLO + SMEFTLO            = " << drEFT.result() << '\n';
  std::cout << "Absolute Delta-r difference   = "
            << absoluteDifference << '\n';
  std::cout << "Relative Delta-r difference   = "
            << relativeDifference << '\n';

  return 0;
}

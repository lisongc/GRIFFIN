/*
  Regression test for smeft::mat_SMEFTLO.

  The test covers the revised local-contact matching for quarks, charged
  leptons, and all neutrino generations; verifies delta S'_4f = 0; checks the
  explicitly linearized Z residue; and confirms that the incomplete nu_e
  full-result path is rejected until its t-channel W module exists.
*/

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>

#include "SMvalG.h"
#include "ff0.h"
#include "smeft/ff0_smeft.h"
#include "smeft/xsc_smeft.h"

using namespace griffin;
using griffin::smeft::EWInputScheme;
using griffin::smeft::Generation;
using griffin::smeft::Input;

namespace {

const Generation one = Generation::first;
const Generation two = Generation::second;
const Generation three = Generation::third;

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

void require(bool condition, const std::string& message)
{
  if(!condition)
    throw std::runtime_error(message);
}

void requireClose(Cplx actual, Cplx expected, const std::string& message)
{
  const double scale = std::max(1.0,std::max(std::abs(actual),std::abs(expected)));
  if(std::abs(actual-expected) > 2e-13*scale)
    throw std::runtime_error(message);
}

template <typename Exception, typename Function>
void requireThrows(Function function, const std::string& message)
{
  bool threw = false;
  try
  {
    function();
  }
  catch(const Exception&)
  {
    threw = true;
  }
  require(threw,message);
}

double projected(double ll, double lr, double rl, double rr,
                 Current initialCurrent, Current finalCurrent,
                 double inverseLambdaSquared)
{
  double combination = 0.0;
  if(initialCurrent == Current::vector && finalCurrent == Current::vector)
    combination = ll+lr+rl+rr;
  else if(initialCurrent == Current::vector && finalCurrent == Current::axial)
    combination = ll-lr+rl-rr;
  else if(initialCurrent == Current::axial && finalCurrent == Current::vector)
    combination = ll+lr-rl-rr;
  else
    combination = ll-lr-rl+rr;
  return 0.25*inverseLambdaSquared*combination;
}

const char* currentName(Current current)
{
  return current == Current::vector ? "V" : "A";
}

void testDownQuarkContact(const inval& smInput)
{
  Input input;
  input.setLambdaGeV(1000.0);
  input.setClq1(one,one,one,one,0.8);
  input.setClq3(one,one,one,one,-0.2);
  input.setCld(one,one,one,one,-0.4);
  input.setCqe(one,one,one,one,0.3);
  input.setCed(one,one,one,one,-0.1);

  FA_SMLO faElectron(Fermion::electron,smInput);
  FA_SMLO faDown(Fermion::d,smInput);
  SW_SMLO swElectron(Fermion::electron,smInput);
  SW_SMLO swDown(Fermion::d,smInput);

  matel sm(Fermion::electron,Fermion::d,Current::vector,Current::vector,
           faElectron,faDown,swElectron,swDown,91.1876*91.1876,0.3,smInput);
  smeft::mat_SMEFTLO eft(
      Fermion::electron,Fermion::d,Current::vector,Current::vector,
      faElectron,faDown,swElectron,swDown,91.1876*91.1876,0.3,smInput,input);

  const Current currents[2] = {Current::vector,Current::axial};
  const double invLambda2 = input.inverseLambdaSquaredGeV();
  for(int initial = 0; initial < 2; ++initial)
  {
    for(int final = 0; final < 2; ++final)
    {
      sm.setform(currents[initial],currents[final]);
      eft.setform(currents[initial],currents[final]);
      const double expected = projected(0.6,-0.4,0.3,-0.1,
                                        currents[initial],currents[final],
                                        invLambda2);
      std::cout << "Down contact " << currentName(currents[initial])
                << currentName(currents[final])
                << " actual=" << eft.coeffS4fD6Tree()
                << " expected=" << expected << '\n';
      requireClose(eft.coeffS4fD6Tree(),expected,
                   "down-quark contact projection has wrong signs");
      requireClose(eft.coeffS()-sm.coeffS(),expected,
                   "down-quark contact was not added to S exactly once");
      requireClose(eft.coeffSp4fD6Tree(),0.0,
                   "local down-quark contact generated nonzero S prime");
      requireClose(eft.coeffSp(),sm.coeffSp(),
                   "SMEFT contact changed the inherited S prime");
      requireClose(eft.resoffZ4fD6Tree(),expected,
                   "exact local background does not equal contact S");
      requireClose(eft.resoffZD6Tree(),expected,
                   "total D6 background does not contain the contact once");
      requireClose(eft.poleD6Tree(),0.0,
                   "contact-only input unexpectedly shifted the pole");
      requireClose(eft.resultD6Tree(),expected,
                   "D6 result is not pole plus exact background");
      requireClose(eft.resoffZ()-sm.resoffZ(),expected,
                   "contact was not added once to the regular amplitude");
      requireClose(eft.result()-sm.result(),expected,
                   "contact-only matrix-element difference is incorrect");
    }
  }

  // The two public constructor styles must be numerically identical.
  smeft::mat_SMEFTLO legacy(
      ELE,DQU,VEC,VEC,faElectron,faDown,swElectron,swDown,
      91.1876*91.1876,0.3,smInput,input);
  eft.setform(Current::vector,Current::vector);
  requireClose(legacy.result(),eft.result(),
               "legacy and enum matrix constructors disagree");

  std::cout << "Down VV contact S [GeV^-2]    = "
            << eft.coeffS4fD6Tree() << '\n';
}

void testAdditionalFlavorAssembly()
{
  const Current currents[2] = {Current::vector,Current::axial};

  Input upInput;
  upInput.setLambdaGeV(1000.0);
  upInput.setClq1(one,one,two,two,0.5);
  upInput.setClq3(one,one,two,two,0.2);
  upInput.setClu(one,one,two,two,0.1);
  // Exercise the old current ordering: Ceq[11ii] = Cqe[ii11].
  upInput.setCeq(one,one,two,two,-0.4);
  upInput.setCeu(one,one,two,two,0.6);
  for(int initial = 0; initial < 2; ++initial)
  {
    for(int final = 0; final < 2; ++final)
    {
      const Cplx actual = smeft::coeffS4fD6Tree(
          Fermion::electron,Fermion::c,currents[initial],currents[final],
          upInput);
      const double expected = projected(0.3,0.1,-0.4,0.6,
                                        currents[initial],currents[final],
                                        upInput.inverseLambdaSquaredGeV());
      requireClose(actual,expected,
                   "up-quark flavor assembly or Ceq alias is incorrect");
    }
  }

  Input tauInput;
  tauInput.setLambdaGeV(1000.0);
  tauInput.setCll(one,one,three,three,0.2);
  tauInput.setCll(one,three,three,one,0.3);
  tauInput.setCle(one,one,three,three,-0.1);
  tauInput.setCle(three,three,one,one,0.4);
  // Enter RR through the Fierz-equivalent alias, while xsc_smeft reads the
  // canonical Cee[1,1,3,3] slot.
  tauInput.setCee(one,three,three,one,0.2);
  for(int initial = 0; initial < 2; ++initial)
  {
    for(int final = 0; final < 2; ++final)
    {
      const Cplx actual = smeft::coeffS4fD6Tree(
          Fermion::electron,Fermion::tau,currents[initial],currents[final],
          tauInput);
      const double expected = projected(0.5,-0.1,0.4,0.2,
                                        currents[initial],currents[final],
                                        tauInput.inverseLambdaSquaredGeV());
      requireClose(actual,expected,
                   "charged-lepton flavor assembly or Cee alias is incorrect");
    }
  }

  requireThrows<std::invalid_argument>([&upInput]() {
    smeft::coeffS4fD6Tree(Fermion::muon,Fermion::c,
                          Current::vector,Current::vector,upInput);
  }, "contact S must reject a non-electron initial state");
  requireThrows<std::domain_error>([&tauInput]() {
    smeft::coeffS4fD6Tree(Fermion::electron,Fermion::electron,
                          Current::vector,Current::vector,tauInput);
  }, "generic contact S must reject Bhabha scattering");
  requireThrows<std::invalid_argument>([&tauInput]() {
    smeft::coeffS4fD6Tree(Fermion::electron,Fermion::tau,
                          Current::scalar,Current::vector,tauInput);
  }, "contact S must reject scalar currents");

  std::cout << "Up/charged-lepton assembly       = PASS\n";
}

void testNeutrinoContacts(const inval& smInput)
{
  FA_SMLO faElectron(Fermion::electron,smInput);
  SW_SMLO swElectron(Fermion::electron,smInput);
  psobsfix neutralSW(0.0,smInput);

  Input muonNeutrinoInput;
  muonNeutrinoInput.setLambdaGeV(1000.0);
  muonNeutrinoInput.setCll(one,one,two,two,0.6);
  muonNeutrinoInput.setCle(two,two,one,one,-0.2);
  // This exchange contraction must not enter a neutral-lepton final state.
  muonNeutrinoInput.setCll(one,two,two,one,9.0);

  FA_SMLO faNuMuon(Fermion::nuMuon,smInput);
  smeft::mat_SMEFTLO nuMuon(
      Fermion::electron,Fermion::nuMuon,Current::vector,Current::vector,
      faElectron,faNuMuon,swElectron,neutralSW,
      91.1876*91.1876,-0.2,smInput,muonNeutrinoInput);

  const double invLambda2 = muonNeutrinoInput.inverseLambdaSquaredGeV();
  requireClose(nuMuon.coeffS4fD6Tree(),
               projected(0.6,0.0,-0.2,0.0,Current::vector,
                         Current::vector,invLambda2),
               "nu_mu VV contact is incorrect");
  nuMuon.setform(Current::vector,Current::axial);
  requireClose(nuMuon.coeffS4fD6Tree(),
               projected(0.6,0.0,-0.2,0.0,Current::vector,
                         Current::axial,invLambda2),
               "nu_mu VA must equal nu_mu VV");
  nuMuon.setform(Current::axial,Current::vector);
  const Cplx nuMuonAV = nuMuon.coeffS4fD6Tree();
  requireClose(nuMuonAV,
               projected(0.6,0.0,-0.2,0.0,Current::axial,
                         Current::vector,invLambda2),
               "nu_mu AV contact is incorrect");
  nuMuon.setform(Current::axial,Current::axial);
  requireClose(nuMuon.coeffS4fD6Tree(),nuMuonAV,
               "nu_mu AA must equal nu_mu AV");
  requireClose(nuMuon.coeffSp4fD6Tree(),0.0,
               "nu_mu local contact generated nonzero S prime");

  Input electronNeutrinoInput;
  electronNeutrinoInput.setLambdaGeV(1000.0);
  electronNeutrinoInput.setCll(one,one,one,one,0.5);
  electronNeutrinoInput.setCle(one,one,one,one,-0.25);

  FA_SMLO faNuElectron(Fermion::nuElectron,smInput);
  smeft::mat_SMEFTLO nuElectron(
      Fermion::electron,Fermion::nuElectron,
      Current::vector,Current::vector,
      faElectron,faNuElectron,swElectron,neutralSW,
      91.1876*91.1876,0.1,smInput,electronNeutrinoInput);

  const double expectedVV = projected(1.0,0.0,-0.25,0.0,
                                      Current::vector,Current::vector,
                                      electronNeutrinoInput
                                          .inverseLambdaSquaredGeV());
  requireClose(nuElectron.coeffS4fD6Tree(),expectedVV,
               "nu_e contact missed the factor two multiplying Cll[1111]");
  requireClose(nuElectron.coeffS(),expectedVV,
               "nu_e S should contain only the local contact in this module");
  requireClose(nuElectron.coeffSp4fD6Tree(),0.0,
               "nu_e local contact generated nonzero S prime");
  requireClose(nuElectron.coeffSp(),0.0,
               "nu_e neutral-current S prime should vanish in this module");
  requireClose(nuElectron.resoffZ4fD6Tree(),expectedVV,
               "nu_e local contact is not available as a regular piece");
  requireThrows<std::domain_error>([&nuElectron]() {
    nuElectron.resoffZD6Tree();
  }, "nu_e background must reject the missing t-channel W contribution");
  requireThrows<std::domain_error>([&nuElectron]() {
    nuElectron.resultD6Tree();
  }, "nu_e D6 result must reject the missing t-channel W contribution");
  requireThrows<std::domain_error>([&nuElectron]() { nuElectron.resoffZ(); },
      "nu_e resoffZ must reject the missing t-channel W background");
  requireThrows<std::domain_error>([&nuElectron]() { nuElectron.result(); },
      "nu_e result must reject the missing t-channel W background");

  std::cout << "Nu_mu AV contact S [GeV^-2]   = " << nuMuonAV << '\n';
  std::cout << "Nu_e VV contact S [GeV^-2]    = " << expectedVV << '\n';
  std::cout << "Nu_e full-result guard         = PASS\n";
}

void testLinearizedResidue(const inval& smInput)
{
  Input input;
  input.setLambdaGeV(1000.0);
  input.setCphiWB(0.012);
  input.setCphiD(-0.018);
  input.setCphil1(one,one,0.007);
  input.setCphil3(one,one,-0.005);
  input.setCphie(one,one,0.004);
  input.setCphiq1(two,two,-0.006);
  input.setCphiq3(two,two,0.009);
  input.setCphiu(two,two,-0.003);
  input.setCphil3(two,two,0.002);
  input.setCll(one,two,two,one,0.001);
  input.setCll(two,one,one,two,0.001);
  // Independent contact contribution for the combined no-double-counting
  // check below.  This four-index coefficient does not enter the Z vertex.
  input.setClq1(one,one,two,two,0.4);

  FA_SMLO faElectron(Fermion::electron,smInput);
  FA_SMLO faCharm(Fermion::c,smInput);
  SW_SMLO swElectron(Fermion::electron,smInput);
  SW_SMLO swCharm(Fermion::c,smInput);

  smeft::mat_SMEFTLO matrix(
      Fermion::electron,Fermion::c,Current::vector,Current::axial,
      faElectron,faCharm,swElectron,swCharm,
      91.1876*91.1876,0.4,smInput,input,EWInputScheme::GmuAlphaMZ);
  matel sm(Fermion::electron,Fermion::c,Current::vector,Current::axial,
           faElectron,faCharm,swElectron,swCharm,
           91.1876*91.1876,0.4,smInput);

  const Cplx ziSM = z0(ELE,VEC,smInput);
  const Cplx zfSM = z0(CQU,AXV,smInput);
  const double dzi = smeft::zD6Tree(
      ELE,VEC,smInput,input,EWInputScheme::GmuAlphaMZ);
  const double dzf = smeft::zD6Tree(
      CQU,AXV,smInput,input,EWInputScheme::GmuAlphaMZ);
  const Cplx expectedD6 = dzi*zfSM+ziSM*dzf;

  requireClose(matrix.matel::coeffR(),ziSM*zfSM,
               "SM form-factor objects do not reconstruct the SM residue");
  requireClose(matrix.coeffRD6Tree(),expectedD6,
               "residue was not linearized with one D6 insertion");
  requireClose(matrix.coeffR(),ziSM*zfSM+expectedD6,
               "cumulative SMEFT residue is incorrect");
  const double mzPole = smInput.get(InputPar::ZMassComplex);
  const double gzPole = smInput.get(InputPar::ZWidthComplex);
  requireClose(matrix.poleD6Tree(),
               expectedD6/Cplx(91.1876*91.1876-
                               mzPole*mzPole,mzPole*gzPole),
               "D6 pole amplitude does not use the complex pole denominator");
  const double expectedBackground = 0.25*0.4
      *input.inverseLambdaSquaredGeV();
  requireClose(matrix.resoffZD6Tree(),expectedBackground,
               "combined test has the wrong exact contact background");
  requireClose(matrix.resultD6Tree(),
               matrix.poleD6Tree()+expectedBackground,
               "D6 result is not exactly pole plus background");
  requireClose(matrix.result()-sm.result(),matrix.resultD6Tree(),
               "cumulative result double-counted or omitted a D6 component");

  std::cout << "Linear residue shift            = " << expectedD6 << '\n';
  std::cout << "Pole/background counting check  = PASS\n";
}

void testZeroLimit(const inval& smInput)
{
  Input input;
  input.setLambdaGeV(1000.0);
  FA_SMLO faElectron(Fermion::electron,smInput);
  FA_SMLO faTau(Fermion::tau,smInput);
  SW_SMLO swElectron(Fermion::electron,smInput);
  SW_SMLO swTau(Fermion::tau,smInput);

  matel sm(Fermion::electron,Fermion::tau,Current::vector,Current::axial,
           faElectron,faTau,swElectron,swTau,90.0*90.0,0.25,smInput);
  smeft::mat_SMEFTLO eft(
      Fermion::electron,Fermion::tau,Current::vector,Current::axial,
      faElectron,faTau,swElectron,swTau,90.0*90.0,0.25,smInput,input);

  requireClose(eft.coeffR(),sm.coeffR(),"zero coefficients changed R");
  requireClose(eft.coeffS(),sm.coeffS(),"zero coefficients changed S");
  requireClose(eft.coeffSp(),sm.coeffSp(),"zero coefficients changed S prime");
  requireClose(eft.poleD6Tree(),0.0,"zero coefficients changed D6 pole");
  requireClose(eft.resoffZD6Tree(),0.0,
               "zero coefficients changed D6 background");
  requireClose(eft.resultD6Tree(),0.0,
               "zero coefficients changed D6 result");
  requireClose(eft.result(),sm.result(),"zero coefficients changed result");
}

void testLinearDifferentialCrossSection(const inval& smInput)
{
  Input plus;
  plus.setLambdaGeV(1000.0);
  plus.setCphiD(0.02);
  plus.setCphiWB(-0.01);
  plus.setCphil3(one,one,0.004);
  plus.setCphil3(two,two,-0.003);
  plus.setCll(one,two,two,one,0.002);
  plus.setCll(two,one,one,two,0.002);
  plus.setCle(one,one,two,two,0.15);

  Input minus = plus;
  minus.clearWilsonCoefficients();
  minus.setCphiD(-0.02);
  minus.setCphiWB(0.01);
  minus.setCphil3(one,one,-0.004);
  minus.setCphil3(two,two,0.003);
  minus.setCll(one,two,two,one,-0.002);
  minus.setCll(two,one,one,two,-0.002);
  minus.setCle(one,one,two,two,-0.15);

  FA_SMLO faElectron(Fermion::electron,smInput);
  FA_SMLO faMuon(Fermion::muon,smInput);
  SW_SMLO swElectron(Fermion::electron,smInput);
  SW_SMLO swMuon(Fermion::muon,smInput);

  const smeft::DifferentialCrossSectionLinear positive =
      smeft::differentialCrossSectionLinear(
          Fermion::electron,Fermion::muon,
          faElectron,faMuon,swElectron,swMuon,
          100.0*100.0,0.3,smInput,plus,
          EWInputScheme::GmuAlphaMZ);
  const smeft::DifferentialCrossSectionLinear negative =
      smeft::differentialCrossSectionLinear(
          Fermion::electron,Fermion::muon,
          faElectron,faMuon,swElectron,swMuon,
          100.0*100.0,0.3,smInput,minus,
          EWInputScheme::GmuAlphaMZ);

  requireClose(positive.sm,negative.sm,
               "SM cross section depends on the SMEFT input");
  requireClose(positive.d6,-negative.d6,
               "D6 cross-section interference is not odd in Wilson coefficients");
  requireClose(positive.total()+negative.total(),2.0*positive.sm,
               "linear cross section retained an O(Lambda^-4) term");
  requireClose(positive.total()-positive.sm,positive.d6,
               "linear cross-section total does not equal SM plus interference");

  std::cout << "Linear cross-section SM [GeV^-2] = " << positive.sm << '\n';
  std::cout << "Linear cross-section D6 [GeV^-2] = " << positive.d6 << '\n';
}

} // namespace

int main()
{
  try
  {
    const SMval directInput = makeDirectMWInput();
    const SMvalGmu smInput(directInput);

    std::cout << std::setprecision(16);
    std::cout << "SMEFT matrix-element test\n";
    std::cout << "=========================\n";
    std::cout << "Input scheme                    = GmuAlphaMZ\n";
    std::cout << "Derived MW pole [GeV]           = "
              << smInput.get(InputPar::WMassComplex) << '\n';

    testDownQuarkContact(smInput);
    testAdditionalFlavorAssembly();
    testNeutrinoContacts(smInput);
    testLinearizedResidue(smInput);
    testZeroLimit(smInput);
    testLinearDifferentialCrossSection(smInput);

    std::cout << "Local four-fermion S prime      = 0\n";
    std::cout << "STATUS                          = PASS\n";
    return 0;
  }
  catch(const std::exception& error)
  {
    std::cerr << "STATUS = FAIL\n";
    std::cerr << error.what() << '\n';
    return 1;
  }
}

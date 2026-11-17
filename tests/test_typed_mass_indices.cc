// Compile in both header orders: mass conventions must not depend on includes.
#ifdef GRIFFIN_SMVAL_FIRST
#include "griffin/sm/SMval.h"
#include "griffin/sm/EWPOZ2.h"
#include "griffin/sm/deltar.h"
#else
#include "griffin/sm/EWPOZ2.h"
#include "griffin/sm/deltar.h"
#include "griffin/sm/SMval.h"
#endif
#include "griffin/sm/SMvalG.h"

#include <cmath>
#include <iostream>
#include <stdexcept>
#include <string>

using namespace griffin;

static_assert(static_cast<int>(InputPar::WMassComplex) == 0, "Pole MW slot changed");
static_assert(static_cast<int>(InputPar::WMassExperiment) == 20, "Running MW slot changed");
static_assert(static_cast<int>(InputPar::ZMassComplex) == 1, "Pole MZ slot changed");
static_assert(static_cast<int>(InputPar::ZMassExperiment) == 21, "Running MZ slot changed");
static_assert(static_cast<int>(InputPar::WWidthComplex) == 17, "Pole GW slot changed");
static_assert(static_cast<int>(InputPar::WWidthExperiment) == 22, "Running GW slot changed");
static_assert(static_cast<int>(InputPar::ZWidthComplex) == 18, "Pole GZ slot changed");
static_assert(static_cast<int>(InputPar::ZWidthExperiment) == 23, "Running GZ slot changed");

#if defined(MW) || defined(MZ) || defined(GamW) || defined(GamZ) \
    || defined(MWc) || defined(MZc) || defined(GWc) || defined(GZc)
#error "Mass headers must not introduce legacy index macros"
#endif

namespace {

void require(bool condition, const char* message)
{
  if(!condition) throw std::runtime_error(message);
}

void close(double actual, double expected, const char* message,
           double relative=2e-12, double absolute=1e-15)
{
  require(std::isfinite(actual) && std::isfinite(expected)
      && std::abs(actual-expected) <= absolute+relative*std::abs(expected), message);
}

void setCouplingsAndMasses(inval& input)
{
  input.set(InputPar::alpha,1.0/137.03599976);
  input.set(InputPar::alphaS,0.1179);
  input.set(InputPar::HiggsMass,125.1);
  input.set(InputPar::topMass,172.5);
  input.set(InputPar::bottomMass,2.87);
  input.set(InputPar::DeltaAlpha,0.059);
  input.set(InputPar::fermiConstant,1.166379e-5);
}

SMval makeInput(double wWidth, double zWidth)
{
  SMval input;
  input.set(InputPar::WMassExperiment,80.377);
  input.set(InputPar::ZMassExperiment,91.1876);
  input.set(InputPar::WWidthExperiment,wWidth);
  input.set(InputPar::ZWidthExperiment,zWidth);
  setCouplingsAndMasses(input);
  return input;
}

void checkConversionsAndPolePrefactors(double wWidth, double zWidth)
{
  const SMval input=makeInput(wWidth,zWidth);
  const double wPole=80.377/std::hypot(1.0,wWidth/80.377);
  const double zPole=91.1876/std::hypot(1.0,zWidth/91.1876);
  close(input.get(InputPar::WMassComplex),wPole,"SMval pole MW conversion");
  close(input.get(InputPar::ZMassComplex),zPole,"SMval pole MZ conversion");
  close(input.get(InputPar::WWidthComplex),wWidth*wPole/80.377,"SMval pole W width");
  close(input.get(InputPar::ZWidthComplex),zWidth*zPole/91.1876,"SMval pole Z width");
  close(input.get(InputPar::WMassExperiment),80.377,"SMval changed running MW",0,0);
  close(input.get(InputPar::ZMassExperiment),91.1876,"SMval changed running MZ",0,0);
  close(input.get(InputPar::WWidthExperiment),wWidth,"SMval changed running W width",0,0);
  close(input.get(InputPar::ZWidthExperiment),zWidth,"SMval changed running Z width",0,0);

  // Independent pole-only input has no running-width slots to fall back to.
  inval poleInput;
  poleInput.set(InputPar::WMassComplex,wPole);
  poleInput.set(InputPar::ZMassComplex,zPole);
  setCouplingsAndMasses(poleInput);
  const double rho=delrho2a2(&poleInput);
  require(std::abs(rho)>1e-8,"The rho correction must exercise the mass prefactors");
  const dr_SMNNLO deltaR(input);
  const SW_SMNNLO weakAngle(Fermion::electron,input);
  close(deltaR.drho2a2(),rho/(1-(zPole/wPole)*(zPole/wPole)),
        "Delta_r inline prefactor must use pole masses in either header order");
  close(weakAngle.drho2a2(),rho*(wPole/zPole)*(wPole/zPole),
        "Weak-angle inline prefactor must use pole masses in either header order");
}

void checkGmuMassSlots()
{
  const SMval running=makeInput(2.085,2.4952);
  const SMvalGmu solved(running);

  inval poleInput;
  poleInput.set(InputPar::ZMassComplex,91.1876/std::hypot(1.0,2.4952/91.1876));
  setCouplingsAndMasses(poleInput);
  const invalGmu poleSolved(poleInput);
  const double wPole=solved.get(InputPar::WMassComplex);
  const double zPole=solved.get(InputPar::ZMassComplex);
  close(wPole,poleSolved.get(InputPar::WMassComplex),"Gmu solve must write pole MW");
  close(zPole,running.get(InputPar::ZMassComplex),"Gmu solve changed pole MZ");
  close(solved.get(InputPar::ZMassExperiment),91.1876,"Gmu solve changed running MZ",0,0);
  close(solved.get(InputPar::ZWidthExperiment),2.4952,"Gmu solve changed running Z width",0,0);
  close(solved.get(InputPar::ZWidthComplex),running.get(InputPar::ZWidthComplex),
        "Gmu solve changed pole Z width");

  const double wPoleWidth=0.3376186*solved.get(InputPar::fermiConstant)
      *wPole*wPole*wPole*(1+0.2122066*solved.get(InputPar::alphaS));
  const double conversion=std::hypot(1.0,wPoleWidth/wPole);
  close(solved.get(InputPar::WWidthComplex),wPoleWidth,"Gmu pole W width");
  close(solved.get(InputPar::WMassExperiment),wPole*conversion,"Gmu running MW");
  close(solved.get(InputPar::WWidthExperiment),wPoleWidth*conversion,"Gmu running W width");
  require(solved.get(InputPar::WMassExperiment)-wPole>0.01,
          "Pole and running MW must remain distinct");

  const dr_SMNNLO deltaR(solved);
  const double rhs=Pi*solved.get(InputPar::alpha)
      /(std::sqrt(2.0)*solved.get(InputPar::fermiConstant))*(1+real(deltaR.result()));
  close(wPole*wPole*(1-wPole*wPole/(zPole*zPole)),rhs,
        "Pole MW must satisfy the Gmu equation",0,2e-4);
}

} // namespace

int main()
{
  try {
    checkConversionsAndPolePrefactors(2.085,2.4952);
    checkConversionsAndPolePrefactors(0.0,0.0);
    checkConversionsAndPolePrefactors(8.0,3.0);
    checkGmuMassSlots();
    std::cout << "Typed mass slots and header order: PASS\n";
    return 0;
  }
  catch(const std::exception& error) {
    std::cerr << "Typed mass slots: FAIL: " << error.what() << '\n';
    return 1;
  }
}

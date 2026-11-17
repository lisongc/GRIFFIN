/*
  File: test_scattering_equivalence.cc
  Purpose: Check legacy SM scattering formulas against the public angular
           helper and the zero-Wilson-coefficient SMEFT calculator.
  Authors: Jackson Wallace and Lisong Chen
  Last updated: 2026-11-16 (migrated into the v1.1.1 SM test group)
*/

// Standalone comparison of the manual calculations in testmatel.cc/testaas.cc
// with the public amplitude-to-cross-section helper and SMEFTCalculator.
// Build/run from the repository root:
//   cmake -S . -B build/scattering-equivalence -DBUILD_TESTING=ON
//   cmake --build build/scattering-equivalence --target test_scattering_equivalence -j 4
//   ./build/scattering-equivalence/tests/test_scattering_equivalence
//
// The calculator always uses mat_SMaas. At physical alpha_s its SM limit
// should match testaas, not the different SM content of testmatel. Both old
// formulas are checked against the helper; at alpha_s=0 the calculator is
// additionally checked against mat_SMNNLO, since the aas extension vanishes.
// This checks API/algebra equivalence using shared SM kernels, not their
// independent physical accuracy. The separate test_near_pole checks the
// corrected upstream v1.1.1 near-pole kernel; exact s=MZ^2 remains unsupported.

#include "griffin/sm/EWPOZ2.h"
#include "griffin/sm/tools.h"
#include "griffin/sm/xscaas.h"
#include "griffin/smeft/SMEFTCalculator.h"
#include "griffin/sm/SMval.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>

using namespace griffin;
using namespace griffin::prediction;

namespace {

typedef std::array<Cplx,4> Amplitudes; // VV, AV, VA, AA
const double toNb = 0.38937966e6;
const double relativeTolerance = 5.0e-12;
const double absoluteTolerance = 1.0e-26;

struct FinalState {
  Fermion fermion;
  const char* name;
  double colors;
};

struct Comparison {
  const char* name;
  unsigned checks = 0;
  double maxRelativeError = 0.0;

  explicit Comparison(const char* label) : name(label) {}

  void check(Cplx actual, Cplx expected, const std::string& point,
             const std::string& quantity)
  {
    const double scale = std::max(std::abs(actual),std::abs(expected));
    const double error = std::abs(actual-expected);
    if(!std::isfinite(real(actual)) || !std::isfinite(imag(actual))
        || !std::isfinite(real(expected)) || !std::isfinite(imag(expected))
        || error > absoluteTolerance+relativeTolerance*scale)
    {
      std::ostringstream message;
      message << std::setprecision(17) << name << ": " << point << ", "
              << quantity << ": old=" << expected << ", new=" << actual
              << ", absolute difference=" << error;
      throw std::runtime_error(message.str());
    }
    ++checks;
    if(scale > 0.0) maxRelativeError = std::max(maxRelativeError,error/scale);
  }

  void report() const
  {
    std::cout << "PASS " << name << ": " << checks
              << " comparisons, max relative difference="
              << std::scientific << std::setprecision(3) << maxRelativeError << '\n';
  }
};

SMval makeInput(double alphaS)
{
  SMval input;
  // Same physical inputs and running-width to pole conversion as the examples.
  input.set(InputPar::alpha,1.0/137.03599976);
  input.set(InputPar::ZMassExperiment,91.1876);
  input.set(InputPar::WMassExperiment,80.377);
  input.set(InputPar::ZWidthExperiment,2.4952);
  input.set(InputPar::WWidthExperiment,2.085);
  input.set(InputPar::HiggsMass,125.1);
  input.set(InputPar::topMass,172.5);
  input.set(InputPar::bottomMass,2.87);
  input.set(InputPar::DeltaAlpha,0.059);
  input.set(InputPar::alphaS,alphaS);
  // The calculator also validates these inputs, omitted by the old examples.
  input.set(InputPar::charmMass,0.652);
  input.set(InputPar::tauMass,1.777);
  input.set(InputPar::strangeMass,0.0);
  input.set(InputPar::upMass,0.0);
  input.set(InputPar::downMass,0.0);
  input.set(InputPar::muonMass,0.0);
  input.set(InputPar::electronMass,0.0);
  input.set(InputPar::fermiConstant,1.166379e-5);
  return input;
}

template<class Matrix>
Amplitudes legacyAmplitudes(Matrix& matrix, double s, double cosine)
{
  matrix.setkinvar(s,cosine);
  Amplitudes result;
  matrix.setform(Current::vector,Current::vector); result[0] = matrix.result();
  matrix.setform(Current::axial,Current::vector); result[1] = matrix.result();
  matrix.setform(Current::vector,Current::axial); result[2] = matrix.result();
  matrix.setform(Current::axial,Current::axial); result[3] = matrix.result();
  return result;
}

double legacyAngularSquare(const Amplitudes& a, double c)
{
  // Explicitly retain the old examples' expression as the reference.
  const Cplx vv=a[0], av=a[1], va=a[2], aa=a[3];
  return real((1+c*c)*(vv*conj(vv)+av*conj(av)+va*conj(va)+aa*conj(aa))
              +4*c*(vv*conj(aa)+va*conj(av)));
}

double helicityMatrixElementSquared(const Amplitudes& a, double s,
                                   double c, double colors)
{
  // Independent algebra: sum the four chiral squares, averaging initial
  // spins and summing final spins/colors. Reduced coefficients have units
  // GeV^-2; the physical |M|^2 is dimensionless and includes Nc*s^2.
  const Cplx ll=a[0]+a[1]+a[2]+a[3], rr=a[0]-a[1]-a[2]+a[3];
  const Cplx lr=a[0]+a[1]-a[2]-a[3], rl=a[0]-a[1]+a[2]-a[3];
  return colors*s*s/4.0*((1+c)*(1+c)*(std::norm(ll)+std::norm(rr))
                        +(1-c)*(1-c)*(std::norm(lr)+std::norm(rl)));
}

void compareHelper(const Amplitudes& old, const FinalState& final,
                   double s, double c, const std::string& context,
                   Comparison& comparison)
{
  const double angular = legacyAngularSquare(old,c);
  const double oldM2 = final.colors*s*s*angular;
  const double oldSigma = final.colors*s/(32*Pi)*angular;
  const double newSigma = differentialCrossSectionFromAmplitudes(
      old,Fermion::electron,final.fermion,s,c);
  comparison.check(helicityMatrixElementSquared(old,s,c,final.colors),oldM2,
                   context,"helicity |M|^2");
  // There is no public |M|^2 accessor: undo the two-body phase-space factor.
  comparison.check(newSigma*(32*Pi*s),oldM2,context,"helper |M|^2");
  comparison.check(newSigma,oldSigma,context,"d sigma/d cos(theta) [GeV^-2]");
  comparison.check(differentialCrossSectionFromAmplitudes(
      old,Fermion::electron,final.fermion,s,c,toNb),oldSigma*toNb,
      context,"d sigma/d cos(theta) [nb]");
  comparison.check(differentialCrossSectionFromAmplitudes(
      old,Fermion::electron,final.fermion,s,c,1000*toNb),oldSigma*(1000*toNb),
      context,"d sigma/d cos(theta) [pb]");
}

void compareCalculator(const Amplitudes& old, const Amplitudes& total,
                       const Amplitudes& sm, const FinitePrediction<double>& sigma,
                       const FinalState& final, double s, double c,
                       const std::string& context, Comparison& comparison)
{
  const char* channels[] = {"VV","AV","VA","AA"};
  for(unsigned k=0;k<old.size();++k)
  {
    comparison.check(total[k],old[k],context,std::string(channels[k])+" amplitude");
    comparison.check(sm[k],old[k],context,std::string(channels[k])+" SM amplitude");
    comparison.check(std::norm(total[k]),std::norm(old[k]),context,
                     std::string(channels[k])+" reduced amplitude square");
  }
  const double angular = legacyAngularSquare(old,c);
  const double oldM2 = final.colors*s*s*angular;
  const double oldSigmaNb = final.colors*s/(32*Pi)*toNb*angular;
  comparison.check(helicityMatrixElementSquared(total,s,c,final.colors),oldM2,
                   context,"calculator helicity |M|^2");
  comparison.check(sigma.total()/toNb*(32*Pi*s),oldM2,
                   context,"calculator cross-section-derived |M|^2");
  comparison.check(sigma.total(),oldSigmaNb,context,"total d sigma/d cos(theta) [nb]");
  comparison.check(sigma.smResult(),oldSigmaNb,context,"SM d sigma/d cos(theta) [nb]");
  comparison.check(sigma.modelCorrection(),0.0,context,"zero-Wilson correction [nb]");
}

unsigned scan(double alphaS, Comparison& nnloHelper, Comparison& aasHelper,
              Comparison& aasCalculator, Comparison& zeroAlphaSCalculator)
{
  const SMval input = makeInput(alphaS);
  const ElectroweakSetup setup(EWInputScheme::AlphaMWMZ,input);
  smeft::Input zero;
  zero.setLambdaGeV(1000.0); // All Wilson coefficients remain zero.
  const SMEFTCalculator calculator(setup,zero);
  const FinalState finals[] = {
    {Fermion::muon,"mu",1}, {Fermion::tau,"tau",1},
    {Fermion::u,"u",3}, {Fermion::d,"d",3}, {Fermion::s,"s",3},
    {Fermion::c,"c",3}, {Fermion::b,"b",3}
  };
  unsigned points = 0;
  for(const FinalState& final : finals)
  {
    const FA_SMNNLO fai(Fermion::electron,input), faf(final.fermion,input);
    const SW_SMNNLO swi(Fermion::electron,input), swf(final.fermion,input);
    // Separate concrete objects preserve the two legacy SM orders.
    mat_SMNNLO nnlo(Fermion::electron,final.fermion,Current::vector,Current::vector,fai,faf,swi,swf,100,.5,input);
    mat_SMaas aas(Fermion::electron,final.fermion,Current::vector,Current::vector,fai,faf,swi,swf,100,.5,input);
    for(double energy=10;energy<=190;energy+=20)
      for(double c : {-0.8,-0.3,0.0,0.5,0.8})
      {
        const double s = energy*energy;
        std::ostringstream label;
        label << "alpha_s=" << alphaS << ", e+e- -> " << final.name
              << " pair, sqrt(s)=" << energy << " GeV, cos(theta)=" << c;
        const std::string context = label.str();
        const Amplitudes oldNNLO = legacyAmplitudes(nnlo,s,c);
        const Amplitudes oldAas = legacyAmplitudes(aas,s,c);
        compareHelper(oldNNLO,final,s,c,context,nnloHelper);
        compareHelper(oldAas,final,s,c,context,aasHelper);

        const ScatteringPoint point(Fermion::electron,final.fermion,s,c,toNb);
        const Amplitudes total = calculator.amplitudes(point);
        const Amplitudes sm = calculator.smAmplitudes(point);
        const auto sigma = calculator.differentialCrossSection(point);
        compareCalculator(oldAas,total,sm,sigma,final,s,c,context,aasCalculator);
        if(alphaS == 0.0)
          compareCalculator(oldNNLO,total,sm,sigma,final,s,c,context,zeroAlphaSCalculator);

        // Show the original examples' d-quark scan alongside the assertions.
        if(alphaS != 0.0 && final.fermion == Fermion::d && c == 0.5)
        {
          const double oldAngular = legacyAngularSquare(oldAas,c);
          std::cout << std::fixed << std::setprecision(0) << std::setw(7) << energy
                    << std::scientific << std::setprecision(9)
                    << "  " << final.colors*s*s*oldAngular
                    << "  " << helicityMatrixElementSquared(total,s,c,final.colors)
                    << "  " << final.colors*s/(32*Pi)*toNb*oldAngular
                    << "  " << sigma.total() << '\n';
        }
        ++points;
      }
  }
  return points;
}

} // namespace

int main()
{
  try
  {
    Comparison nnloHelper("mat_SMNNLO manual vs helper");
    Comparison aasHelper("mat_SMaas manual vs helper");
    Comparison aasCalculator("mat_SMaas manual vs zero-Wilson calculator");
    Comparison zeroAlphaSCalculator("mat_SMNNLO vs calculator at alpha_s=0");
    std::cout << "Legacy/new scattering equivalence; relative tolerance="
              << relativeTolerance << ", absolute tolerance=" << absoluteTolerance
              << "\nPhysical |M|^2 is spin averaged and color summed."
              << "\nSample: e+e- -> d dbar, cos(theta)=0.5, alpha_s=0.1179, mat_SMaas"
              << "\nsqrt(s)   old |M|^2       new |M|^2       old dsigma [nb]  new dsigma [nb]\n";
    unsigned points = scan(0.1179,nnloHelper,aasHelper,aasCalculator,zeroAlphaSCalculator);
    points += scan(0.0,nnloHelper,aasHelper,aasCalculator,zeroAlphaSCalculator);
    std::cout << '\n';
    nnloHelper.report();
    aasHelper.report();
    aasCalculator.report();
    zeroAlphaSCalculator.report();
    std::cout << "PASS: " << points
              << " points (7 final states x 10 energies x 5 angles x 2 alpha_s values).\n";
    return 0;
  }
  catch(const std::exception& error)
  {
    std::cerr << "FAIL: scattering equivalence: " << error.what() << '\n';
    return 1;
  }
}

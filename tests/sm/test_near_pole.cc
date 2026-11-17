/*
  File: test_near_pole.cc
  Purpose: Regress the ported v1.1.1 mixed-background correction against
           independently archived upstream numbers near the real Z pole.
  Authors: Jackson Wallace and Lisong Chen
  Last updated: 2026-11-16
*/

#include "griffin/SMEFT.h"
#include "griffin/sm/ff0.h"
#include "griffin/sm/tools.h"
#include "griffin/sm/xscaas.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <limits>
#include <stdexcept>
#include <string>

using namespace griffin;
using namespace griffin::prediction;

namespace {

void require(bool condition, const std::string& message)
{
  if(!condition) throw std::runtime_error(message);
}

void close(double actual, double expected, double absolute,
           const std::string& message)
{
  // Near-pole subtractions lose digits; these absolute tolerances accommodate
  // floating-point/compiler differences without hiding the old large defect.
  require(std::isfinite(actual) && std::fabs(actual-expected)<=absolute, message);
}

SMval makeInput()
{
  SMval input;
  input.set(InputPar::alpha,1.0/137.03599976);
  input.set(InputPar::alphaS,0.1179);
  input.set(InputPar::ZMassExperiment,91.1876);
  input.set(InputPar::ZWidthExperiment,2.4952);
  input.set(InputPar::WMassExperiment,80.377);
  input.set(InputPar::WWidthExperiment,2.085);
  input.set(InputPar::HiggsMass,125.1);
  input.set(InputPar::topMass,172.5);
  input.set(InputPar::bottomMass,2.87);
  input.set(InputPar::charmMass,0.652);
  input.set(InputPar::tauMass,1.777);
  for(InputPar key : {InputPar::electronMass,InputPar::muonMass,
      InputPar::upMass,InputPar::downMass,InputPar::strangeMass})
    input.set(key,0.0);
  input.set(InputPar::DeltaAlpha,0.059);
  input.set(InputPar::fermiConstant,1.166379e-5);
  return input;
}

struct Result {
  double smNb, mixedNb;
  Cplx smVV, mixedVV;
};

Result evaluate(const SMval& input, Fermion final, double energy)
{
  const double s=energy*energy;
  FA_SMNNLO fae(Fermion::electron,input), faf(final,input);
  SW_SMNNLO swe(Fermion::electron,input), swf(final,input);
  std::array<Cplx,4> a,b;
  Result result;
  for(unsigned k=0;k<4;++k) {
    const Current ci=k%2 ? Current::axial : Current::vector;
    const Current cf=k/2 ? Current::axial : Current::vector;
    const mat_SMNNLO sm(Fermion::electron,final,ci,cf,fae,faf,swe,swf,s,.5,input);
    const mat_SMaas mixed(Fermion::electron,final,ci,cf,fae,faf,swe,swf,s,.5,input);
    a[k]=sm.result(); b[k]=mixed.result();
    if(k==0) {
      result.smVV=sm.resoffZ();
      result.mixedVV=mixed.resoffZ();
    }
  }
  result.smNb=differentialCrossSectionFromAmplitudes(
      a,Fermion::electron,final,s,.5,.38937966e6);
  result.mixedNb=differentialCrossSectionFromAmplitudes(
      b,Fermion::electron,final,s,.5,.38937966e6);
  return result;
}

// Independent source: exact upstream v1.1.1 tag
// a18f53830283a78155750a23a3594ada8de4c074, evaluated 2026-11-14 by
// diagnostics/colleague_near_pole_20261114/compare_near_pole.cc (outside this
// package); v111_email_grid.dat and v111_muon_fine.dat. These constants are
// copied, not generated using the library under test. The 20 down-quark rows
// reproduce the collaborator's printed table. All cross sections are nb,
// cos(theta)=0.5, with final-state color included.
void testUpstreamNumbers(const SMval& input)
{
  // E=91.14+i*0.001 GeV; sigma_SM, sigma_mixed, Re/Im Bvv_SM, Re/Im Bvv_mixed.
  const double down[][6]={
    {4.59339252818531,4.59367149628642,3.9294309395629e-6,-6.61804141968155e-8,3.92985512863177e-6,-6.73512511329924e-8},
    {4.59366879956021,4.59394763005024,3.92934395830531e-6,-6.61713343775216e-8,3.92977157103219e-6,-6.73388948604739e-8},
    {4.59393918024557,4.59421782892202,3.92925695246632e-6,-6.61622543051614e-8,3.92968858445895e-6,-6.7325972629765e-8},
    {4.59420366879104,4.59448207889886,3.9291699200014e-6,-6.61531741702088e-8,3.92960633735595e-6,-6.73123227021277e-8},
    {4.594462263635,4.59474036057025,3.92908285898121e-6,-6.61440938962412e-8,3.92952507021365e-6,-6.72977147518278e-8},
    {4.59471496350413,4.59499264639568,3.92899576682251e-6,-6.61350134544187e-8,3.92944513720922e-6,-6.72818099004781e-8},
    {4.59496176686062,4.59523889452502,3.92890864055124e-6,-6.61259329101798e-8,3.92936708289419e-6,-6.72640881760007e-8},
    {4.59520267252514,4.59547903858674,3.92882147387055e-6,-6.61168522151487e-8,3.92929178644374e-6,-6.72437085123626e-8},
    {4.59543767907959,4.59571296438532,3.92873426492194e-6,-6.61077715465925e-8,3.92922077772837e-6,-6.72192162570259e-8},
    {4.59566678546945,4.59594045741771,3.92864700148555e-6,-6.60986906929291e-8,3.92915694368338e-6,-6.71878575929244e-8},
    {4.59588999012128,4.59616105735246,3.92855967665299e-6,-6.60896097187295e-8,3.929106508325e-6,-6.71437144784187e-8},
    {4.59610729236584,4.59637356502965,3.92847226695464e-6,-6.60805279713446e-8,3.92908572647811e-6,-6.70713252046444e-8},
    {4.59631869126341,4.5965735432003,3.92838474859969e-6,-6.60714494506132e-8,3.92915482963839e-6,-6.69134635516293e-8},
    {4.59652416537801,4.59671961310501,3.92829704037254e-6,-6.60623353112384e-8,3.92987547237252e-6,-6.6136556371014e-8},
    {4.59672379835634,4.59708901811192,3.92820880247702e-6,-6.60537117204646e-8,3.92748285483023e-6,-6.83165982883364e-8},
    {4.59691754175293,4.59722920102186,3.92812115608322e-6,-6.6045476082978e-8,3.92812414907084e-6,-6.76159932080889e-8},
    {4.59710537133546,4.59740607533032,3.9280336325622e-6,-6.60372314748568e-8,3.92818690340097e-6,-6.74649906558446e-8},
    {4.59728729270303,4.59758332859695,3.92794621084846e-6,-6.60289874373376e-8,3.92816436203063e-6,-6.73950988714738e-8},
    {4.59746330432069,4.59775678556286,3.92785886346103e-6,-6.60207433193842e-8,3.92811318500584e-6,-6.73524756272669e-8},
    {4.59763340530264,4.5979252964015,3.92777157224827e-6,-6.60124986891588e-8,3.92804895917558e-6,-6.73222989268895e-8}
  };
  for(unsigned i=0;i<20;++i) {
    const Result value=evaluate(input,Fermion::d,91.14+i*.001);
    const std::string label="upstream down row "+std::to_string(i);
    close(value.smNb,down[i][0],2e-8,label+" SM cross section");
    close(value.mixedNb,down[i][1],2e-8,label+" mixed cross section");
    close(real(value.smVV),down[i][2],5e-13,label+" Re SM VV");
    close(imag(value.smVV),down[i][3],5e-14,label+" Im SM VV");
    close(real(value.mixedVV),down[i][4],5e-13,label+" Re mixed VV");
    close(imag(value.mixedVV),down[i][5],5e-14,label+" Im mixed VV");
  }
  // offset from exact pole in MeV, sigma_SM, sigma_mixed.
  const double muon[][3]={
    {-10,0.95382959408611,0.953858092606414},
    {-1,0.954912072495395,0.954936794293772},
    {-.5,0.954969316585838,0.954989643534945},
    {-.1,0.955006776829663,0.954992154586238},
    {.1,0.955029440547155,0.955102955855666},
    {.5,0.955082907554534,0.95512087817669},
    {1,0.95513925508985,0.955172829841226},
    {10,0.956101291670721,0.95613108593217}
  };
  const double mz=input.get(InputPar::ZMassComplex);
  for(const auto& row : muon) {
    const Result value=evaluate(input,Fermion::muon,mz+row[0]*.001);
    const std::string label="upstream muon offset "+std::to_string(row[0]);
    close(value.smNb,row[1],2e-8,label+" SM cross section");
    close(value.mixedNb,row[2],2e-8,label+" mixed cross section");
  }
}

void testQuadraticBackground(const SMval& input)
{
  const double mz=input.get(InputPar::ZMassComplex), Z=mz*mz;
  const double R0=az0(Fermion::electron,input)*az0(Fermion::muon,input);
  FA_SMNNLO fae(Fermion::electron,input), fam(Fermion::muon,input);
  SW_SMNNLO swe(Fermion::electron,input), swm(Fermion::muon,input);
  Cplx first;
  for(double h : {0.03,0.01,0.003}) {
    const mat_SMaas plus(Fermion::electron,Fermion::muon,Current::axial,
        Current::axial,fae,fam,swe,swm,Z+h,.5,input);
    const mat_SMaas minus(Fermion::electron,Fermion::muon,Current::axial,
        Current::axial,fae,fam,swe,swm,Z-h,.5,input);
    const Cplx bp=plus.resoffZaas(), bm=minus.resoffZaas();
    if(h==.03) first=bp;
    // AA added background is -R0*H(s)/Delta^2, so the upstream quadratic
    // interpolation makes it constant throughout this window, on both sides.
    // H uses a stored MZ^2, while resoffZaas contains s-mz*mz. A compiler
    // may fuse the latter multiply/subtraction: their Delta values can differ
    // by O(epsilon*MZ^2), amplified by division by |h| and by the square.
    // Bound those few operations rather than demand an unconditioned 1e-20.
    // Disabling FMA contraction in xscaas.cc was checked separately and does
    // recover that tighter equality; production arithmetic is unchanged.
    const double roundoff=8*std::numeric_limits<double>::epsilon()
        *(1+Z/h)*std::abs(first);
    if(!(std::abs(bp-first)<roundoff && std::abs(bm-first)<roundoff)) {
      std::cerr<<std::setprecision(17)<<"AA constancy at h="<<h<<": first="<<first
               <<", plus="<<bp<<", minus="<<bm<<", delta="<<bp-first
               <<", minusDelta="<<bm-first<<'\n';
      throw std::runtime_error("mixed AA numerator is not quadratic near the pole");
    }
    const Cplx Q=-h/(2*R0)*(plus.resoffZ()-minus.resoffZ());
    require(std::isfinite(real(Q)) && std::isfinite(imag(Q)),"nonfinite full-background Q");
    // Old imaginary Q was 6.7379973e-4. The corrected contribution vanishes;
    // the inherited SM background leaves only a tiny finite-difference value.
    require(std::fabs(imag(Q))<1e-10,"old imaginary single-pole diagnostic returned");
  }
}

void testExactPoleGuard(const SMval& input)
{
  const ElectroweakSetup setup(EWInputScheme::AlphaMWMZ,input);
  smeft::Input coefficients;
  coefficients.setLambdaGeV(246.0);
  const SMEFTCalculator calculation(setup,coefficients);
  const double mz=calculation.resolvedInput().get(InputPar::ZMassComplex);
  bool rejected=false;
  try {
    calculation.differentialCrossSection(
        ScatteringPoint(Fermion::electron,Fermion::muon,mz*mz,.5));
  } catch(const std::domain_error& error) {
    rejected=std::string(error.what()).find("exact-real-Z-pole")!=std::string::npos;
  }
  require(rejected,"exact real pole must remain explicitly unsupported");
}

} // namespace

int main()
{
  try {
    const SMval input=makeInput();
    testUpstreamNumbers(input);
    testQuadraticBackground(input);
    testExactPoleGuard(input);
    std::cout<<"Near-pole regression passed: 120 down-quark values, 16 muon "
             <<"cross sections, quadratic AA background, exact-pole guard.\n";
    return 0;
  } catch(const std::exception& error) {
    std::cerr<<error.what()<<'\n';
    return 1;
  }
}

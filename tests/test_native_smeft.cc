/*
  File: test_native_smeft.cc
  Purpose: Test finite SMEFT mass matching, direct shifts, observables, amplitude
           assembly, full cross-section squares, ownership and domain validation.
  Authors: Jackson Wallace and Lisong Chen
  Last updated: 2026-11-16
*/

#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <utility>

#include "griffin/sm/EWPOZ2.h"
#include "griffin/sm/SMvalG.h"
#include "griffin/sm/ff.h"
#include "griffin/sm/ff0.h"
#include "griffin/smeft/SMEFTCalculator.h"
#include "griffin/smeft/ff0_smeft.h"
#include "griffin/sm/xscaas.h"

using namespace griffin;
using namespace griffin::prediction;

namespace {

const smeft::Generation one = smeft::Generation::first;
const smeft::Generation two = smeft::Generation::second;
const double pbPerGeV2 = 0.38937966e9;

void require(bool condition, const std::string& message)
{
  if(!condition) throw std::runtime_error(message);
}

void close(Cplx actual, Cplx expected, const std::string& message,
           double relative = 3.0e-12, double absolute = 1.0e-15)
{
  if(!std::isfinite(real(actual)) || !std::isfinite(imag(actual))
      || std::abs(actual-expected)
         > absolute+relative*std::max(std::abs(actual),std::abs(expected)))
    throw std::runtime_error(message);
}

template<typename Callable>
void throws(const Callable& call, const std::string& message)
{
  try { call(); }
  catch(const std::exception&) { return; }
  throw std::runtime_error(message);
}

template<typename Callable>
void throwsInvalidArgument(const Callable& call, const std::string& message)
{
  try { call(); }
  catch(const std::invalid_argument&) { return; }
  throw std::runtime_error(message);
}

SMval makeSM()
{
  SMval input;
  input.set(MZ,91.1876); input.set(MW,80.377);
  input.set(GamZ,2.4952); input.set(GamW,2.085);
  input.set(al,1.0/137.03599976); input.set(als,0.1179);
  input.set(MH,125.1); input.set(MT,172.5); input.set(MB,2.87);
  input.set(MC,0.652); input.set(ML,1.777);
  input.set(MS,0.0); input.set(MU,0.0); input.set(MD,0.0);
  input.set(MM,0.0); input.set(ME,0.0);
  input.set(Delal,0.059); input.set(Gmu,1.166379e-5);
  return input;
}

smeft::Input zero(double lambda = 246.0)
{
  smeft::Input result;
  result.setLambdaGeV(lambda);
  return result;
}

smeft::Input benchmark(double factor = 1.0, double lambda = 246.0)
{
  smeft::Input result = zero(lambda);
  result.setCphiD(0.0003*factor);
  result.setCphiWB(0.0002*factor);
  result.setCphil1(one,one,0.0001*factor);
  result.setCphil3(two,two,-0.00015*factor);
  result.setCphie(two,two,0.00012*factor);
  result.setCll(one,one,two,two,0.0004*factor);
  result.setCll(one,two,two,one,0.0001*factor);
  result.setCle(one,one,two,two,-0.0002*factor);
  result.setCle(two,two,one,one,0.0003*factor);
  result.setCee(one,one,two,two,0.0002*factor);
  return result;
}

ScatteringPoint point(double energy = 90.0, Fermion final = Fermion::muon)
{
  return ScatteringPoint(Fermion::electron,final,energy*energy,0.5,pbPerGeV2);
}

SMEFTCalculator::Amplitudes independentSM(
    const inval& input, const ScatteringPoint& where)
{
  const FA_SMNNLO fai(where.initial,input), faf(where.final,input);
  const SW_SMNNLO swi(where.initial,input), swf(where.final,input);
  const Current initial[]={Current::vector,Current::axial,Current::vector,Current::axial};
  const Current final[]={Current::vector,Current::vector,Current::axial,Current::axial};
  SMEFTCalculator::Amplitudes result;
  for(std::size_t k=0; k<result.size(); ++k)
  {
    const mat_SMaas matrix(where.initial,where.final,
        initial[k],final[k],
        fai,faf,swi,swf,where.s,where.cosTheta,input);
    result[k] = matrix.result();
  }
  return result;
}

double independentSquare(const SMEFTCalculator::Amplitudes& a,
                          const ScatteringPoint& where)
{
  const double c = where.cosTheta;
  const double color = (where.final == Fermion::d || where.final == Fermion::u
      || where.final == Fermion::s || where.final == Fermion::c
      || where.final == Fermion::b) ? 3.0 : 1.0;
  double even = 0.0;
  for(Cplx value : a) even += std::norm(value);
  return color*where.s/(32.0*Pi)*where.unitConversion
    *((1.0+c*c)*even+4.0*c*real(a[0]*conj(a[3])+a[2]*conj(a[1])));
}

void samePrediction(const SMEFTCalculator& first,
                    const SMEFTCalculator& second,
                    const std::string& label)
{
  for(PseudoObservable observable : {PseudoObservable::SW,PseudoObservable::FA,
                                     PseudoObservable::FV,PseudoObservable::WMass})
    close(first.pseudoObservable(observable,Fermion::electron).total(),
          second.pseudoObservable(observable,Fermion::electron).total(),label);
  const auto a = first.amplitudes(point()), b = second.amplitudes(point());
  for(std::size_t k=0; k<a.size(); ++k) close(a[k],b[k],label);
  close(first.differentialCrossSection(point()).total(),
        second.differentialCrossSection(point()).total(),label);
}

void testZeroAndSquare(const ElectroweakSetup& context)
{
  const SMEFTCalculator prediction(context,zero());
  const inval& input = prediction.resolvedInput();
  close(input.get(InputPar::WMassComplex),
        context.input().get(InputPar::WMassComplex),"zero SMEFT changed SM MW",
        2.0e-11,1.0e-10);
  for(Fermion fermion : {Fermion::electron,Fermion::muon,Fermion::tau,
                         Fermion::u,Fermion::b})
  {
    close(prediction.pseudoObservable(PseudoObservable::SW,fermion).total(),
          SW_SMNNLO(fermion,input).result(),"zero SW is not highest SM");
    close(prediction.pseudoObservable(PseudoObservable::FA,fermion).total(),
          FA_SMNNLO(fermion,input).result(),"zero FA is not highest SM");
    // Legacy FV constructors allocate their own FA/SW helpers. This small
    // comparison deliberately tests the existing public highest-SM result.
    close(prediction.pseudoObservable(PseudoObservable::FV,fermion).total(),
          FV_SMNNLO(fermion,input).result(),"zero FV is not highest SM");
    close(prediction.pseudoObservable(PseudoObservable::FA,fermion).modelCorrection(),
          0.0,"zero FA has a model correction",0.0,1.0e-15);
  }
  for(Fermion final : {Fermion::muon,Fermion::b})
    for(double energy : {50.0,90.0,95.0,120.0})
    {
      const ScatteringPoint where = point(energy,final);
      const auto expected = independentSM(input,where);
      const auto actual = prediction.amplitudes(where);
      const auto contacts = prediction.contactAmplitudes(where);
      for(std::size_t k=0; k<actual.size(); ++k)
      {
        close(actual[k],expected[k],"zero amplitude differs from full SM");
        close(contacts[k],0.0,"zero input produced contact amplitude",0.0,0.0);
      }
      const auto result = prediction.differentialCrossSection(where);
      close(result.total(),independentSquare(expected,where),"full SM square mismatch");
      close(result.smResult(),result.total(),"zero model differs from SM result");
    }

  const SMEFTCalculator nonzero(context,benchmark());
  const Current initial[]={Current::vector,Current::axial,Current::vector,Current::axial};
  const Current final[]={Current::vector,Current::vector,Current::axial,Current::axial};
  for(double energy : {50.0,90.0,120.0})
  {
    const ScatteringPoint where = point(energy);
    const auto a = nonzero.amplitudes(where);
    close(nonzero.differentialCrossSection(where).total(),
          independentSquare(a,where),"nonzero prediction is not the full square");
    for(std::size_t k=0; k<a.size(); ++k)
    {
      const CurrentChannel channel(initial[k],final[k]);
      close(nonzero.amplitude(where,channel).total(),a[k],"single-channel ordering");
    }
  }
}

void testFVAndLambda(const ElectroweakSetup& context)
{
  const SMEFTCalculator prediction(context,benchmark());
  const SMEFTCalculator scaled(context,benchmark(4.0,492.0));
  samePrediction(prediction,scaled,"C/Lambda^2 rescaling changed finite result");
  const inval& input = prediction.resolvedInput();
  for(Fermion fermion : {Fermion::electron,Fermion::muon,Fermion::b})
  {
    const int type = index(fermion);
    const Cplx fa = prediction.pseudoObservable(PseudoObservable::FA,fermion).total();
    const Cplx sw = prediction.pseudoObservable(PseudoObservable::SW,fermion).total();
    const Cplx fv = prediction.pseudoObservable(PseudoObservable::FV,fermion).total();
    const double a0 = az0(type,input), v0 = vz0(type,input);
    const double imaginaryRatio =
      (a0*(iz1f(type,VEC,input)+iz1b(type,VEC,input))
       -v0*(iz1f(type,AXV,input)+iz1b(type,AXV,input)))/(a0*a0);
    const double q = 1.0-4.0*std::fabs(Qf[type])*real(sw);
    close(fv,fa*(q*q+imaginaryRatio*imaginaryRatio),"FV lost FA/SW identity");
    require(std::abs(fv-fa*q*q)>1.0e-15,"FV test has no absorptive contribution");
    const double dv = smeft::directZVectorShift(type,input,benchmark());
    require(std::abs(fv-(v0*v0+2.0*v0*dv))>1.0e-8,
            "FV reverted to independent tree direct formula");
  }
}

// Independently evaluate the native finite fixed-point equation at its
// returned point. This checks the equation itself, not a derivative or a
// first-order approximation to its solution.
void checkFiniteEquation(const SMEFTCalculator& prediction,
                         const smeft::WMassSolverOptions& options)
{
  const inval& input=prediction.resolvedInput();
  const smeft::Input& coefficients=prediction.wilsonInput();
  const double mw=input.get(InputPar::WMassComplex);
  const double mz=input.get(InputPar::ZMassComplex);
  const double gmu=input.get(InputPar::fermiConstant);
  const double alpha=input.get(InputPar::alpha);
  const double k=Pi*alpha/(std::sqrt(2.0)*gmu*mz*mz);
  const double root=std::sqrt(1.0-4.0*k);
  const double cw2=(1.0+root)/2.0, sw2=(1.0-root)/2.0;
  const double muonDecay=coefficients.Cphil3(one,one)
      +coefficients.Cphil3(two,two)
      -.5*coefficients.Cll(one,two,two,one);
  // This is a unique-operator strength. The reversed spelling is an alias,
  // so adding its getter here would count the same operator twice.
  const double matching=coefficients.inverseLambdaSquaredGeV()
      /(std::sqrt(2.0)*gmu)
      *(2.0*std::sqrt(cw2/sw2)*coefficients.CphiWB()
         +cw2/(2.0*sw2)*coefficients.CphiD()+muonDecay);
  const double additive=-k/root*matching;
  close(prediction.massSolution().additiveMassSquaredShiftOverMZ2(),additive,
        "finite solver inserted the wrong muon-decay matching term",2e-12,1e-17);
  const double smDeltaR=real(dr_SMNNLO(input).result());
  const double next=mz*mz*(.5+.5*std::sqrt(1.0-4.0*k*(1.0+smDeltaR))+additive);
  const double residual=std::fabs(next-mw*mw);
  // Independent algebra/rounding of the square root and MW*MW can differ by
  // a few double-precision ulps from the stored iteration's arithmetic.
  require(residual<=options.massSquaredToleranceGeV2+5e-11,
          "returned MW does not solve the independently evaluated finite equation");
  close(residual,prediction.massSolution().finalResidualGeV2(),
        "reported finite solver residual is inconsistent",0.0,5e-11);
}

void testFiniteGmu(const ElectroweakSetup& context)
{
  smeft::WMassSolverOptions tight;
  tight.massSquaredToleranceGeV2 = 1.0e-10;
  smeft::Input plus = zero(), minus = zero();
  plus.setCphiWB(0.00060516); minus.setCphiWB(-0.00060516);
  const SMEFTCalculator p(context,plus,tight), m(context,minus,tight);
  const SMEFTCalculator z(context,zero(),tight);
  const double mwp = p.resolvedInput().get(InputPar::WMassComplex);
  const double mwm = m.resolvedInput().get(InputPar::WMassComplex);
  const double mw0 = z.resolvedInput().get(InputPar::WMassComplex);
  require(mwp<mw0 && mwm>mw0,"finite Gmu MW has incorrect shift direction");
  require(std::fabs(mwp+mwm-2.0*mw0)>1.0e-7,
          "finite Gmu MW was accidentally linearized");
  const SMEFTCalculator mixed(context,benchmark(),tight);
  checkFiniteEquation(p,tight); checkFiniteEquation(m,tight);
  checkFiniteEquation(z,tight); checkFiniteEquation(mixed,tight);
  require(p.massSolution().iterations()>0,"Gmu input was not solved");
  require(p.massSolution().finalResidualGeV2()<=tight.massSquaredToleranceGeV2,
          "Gmu convergence residual exceeds tolerance");
  const double gw = p.resolvedInput().get(InputPar::WWidthComplex);
  const double expectedWidth = 0.3376186*p.resolvedInput().get(InputPar::fermiConstant)
      *powint(mwp,3)*(1.0+0.2122066*p.resolvedInput().get(InputPar::alphaS));
  close(gw,expectedWidth,"finite mass did not update W width");
  close(p.resolvedInput().get(InputPar::WMassExperiment),
        mwp*std::sqrt(1.0+sqr(gw/mwp)),"running/pole W conversion changed");
}

void testContactOnce(const ElectroweakSetup& context)
{
  smeft::Input plus = zero(), minus = zero();
  plus.setCee(one,one,two,two,0.002);
  minus.setCee(one,one,two,two,-0.002);
  const SMEFTCalculator p(context,plus), m(context,minus), z(context,zero());
  const ScatteringPoint where = point(50.0);
  const auto a = p.amplitudes(where), b = m.amplitudes(where);
  const auto sm = z.amplitudes(where), contact = p.contactAmplitudes(where);
  const double strength = 0.002/(4.0*246.0*246.0);
  const double signs[4] = {1.0,-1.0,-1.0,1.0};
  for(std::size_t k=0; k<a.size(); ++k)
  {
    close(contact[k],signs[k]*strength,"canonical RR contact projection changed");
    close(a[k]-sm[k],contact[k],"contact not added exactly once");
    close(a[k]+b[k],2.0*sm[k],"contact amplitude is not linear");
  }
  const double quadratic = independentSquare(contact,where);
  const double even = p.differentialCrossSection(where).modelCorrection()
      +m.differentialCrossSection(where).modelCorrection();
  require(quadratic>0.0,"contact-square test is degenerate");
  close(even,2.0*quadratic,"full square lost its contact quadratic term",
        1.0e-6,1.0e-12);
  require(even>0.0,"contact-only full prediction was linearized");
}

void sameContactPredictions(const SMEFTCalculator& actual,
                            const SMEFTCalculator& expected,
                            const std::string& label)
{
  samePrediction(actual,expected,label);
  for(double energy:{50.0,90.0,95.0,120.0})
  {
    const ScatteringPoint where=point(energy);
    const auto a=actual.amplitudes(where), b=expected.amplitudes(where);
    const auto ac=actual.contactAmplitudes(where), bc=expected.contactAmplitudes(where);
    for(std::size_t k=0;k<a.size();++k)
    {
      close(a[k],b[k],label+" amplitude");
      close(ac[k],bc[k],label+" contact");
    }
    const auto sa=actual.differentialCrossSection(where);
    const auto sb=expected.differentialCrossSection(where);
    close(sa.smResult(),sb.smResult(),label+" SM result");
    close(sa.total(),sb.total(),label+" full cross section");
    close(sa.modelCorrection(),sb.modelCorrection(),label+" correction");
  }
}

void testEquivalentContactLabels(const ElectroweakSetup& context)
{
  const smeft::Generation rr[4][4]={{one,one,two,two},{two,two,one,one},
      {one,two,two,one},{two,one,one,two}};
  const smeft::Generation ll[2][2][4]={
      {{one,one,two,two},{two,two,one,one}},
      {{one,two,two,one},{two,one,one,two}}
  };
  smeft::WMassSolverOptions tight;
  tight.massSquaredToleranceGeV2=1e-10;
  const double c=.00060516;
  const SMEFTCalculator sm(context,zero(),tight);

  // A lone RR coefficient has four equivalent setter spellings, not four
  // independently added strengths. Repeating every spelling remains one C.
  smeft::Input rrCanonical=zero(), rrRepeated=zero();
  rrCanonical.setCee(one,one,two,two,c);
  const SMEFTCalculator expectedRR(context,rrCanonical,tight);
  for(const auto& label:rr)
  {
    smeft::Input input=zero();
    input.setCee(label[0],label[1],label[2],label[3],c);
    rrRepeated.setCee(label[0],label[1],label[2],label[3],c);
    sameContactPredictions(SMEFTCalculator(context,input,tight),expectedRR,
                           "single equivalent Cee spelling");
  }
  sameContactPredictions(SMEFTCalculator(context,rrRepeated,tight),expectedRR,
                         "repeated Cee spellings must not accumulate");

  for(int contraction=0;contraction<2;++contraction)
  {
    const auto& canonical=ll[contraction][0];
    smeft::Input input=zero(), repeated=zero();
    input.setCll(canonical[0],canonical[1],canonical[2],canonical[3],c);
    const SMEFTCalculator expected(context,input,tight);
    for(const auto& label:ll[contraction])
    {
      smeft::Input alias=zero();
      alias.setCll(label[0],label[1],label[2],label[3],c);
      repeated.setCll(label[0],label[1],label[2],label[3],c);
      sameContactPredictions(SMEFTCalculator(context,alias,tight),expected,
                             "single equivalent Cll spelling");
    }
    sameContactPredictions(SMEFTCalculator(context,repeated,tight),expected,
                           "repeated Cll pair must not accumulate");
    const double mw=expected.resolvedInput().get(InputPar::WMassComplex);
    const double mwSM=sm.resolvedInput().get(InputPar::WMassComplex);
    if(context.inputScheme()==EWInputScheme::GmuAlphaMZ && contraction==1)
    {
      // Independent one-coefficient normalization, with no sum of aliased
      // getters: delta_mu = -c/2, exactly once.
      const inval& resolved=expected.resolvedInput();
      const double mz=resolved.get(InputPar::ZMassComplex);
      const double gmu=resolved.get(InputPar::fermiConstant);
      const double k=Pi*resolved.get(InputPar::alpha)/(std::sqrt(2.0)*gmu*mz*mz);
      const double epsilon=1.0/(std::sqrt(2.0)*gmu*246.0*246.0);
      const double expectedAdditive=k/std::sqrt(1.0-4.0*k)*epsilon*c/2.0;
      close(expected.massSolution().additiveMassSquaredShiftOverMZ2(),expectedAdditive,
            "canonical Cll muon-decay strength was not counted once",2e-12,1e-17);
      require(mw>mwSM,"positive canonical Cll1221 has no positive finite MW shift");
      checkFiniteEquation(expected,tight);
    }
    else
    {
      close(mw,mwSM,"Cll1122 or direct-MW contact unexpectedly shifted MW",0.0,0.0);
      close(expected.massSolution().additiveMassSquaredShiftOverMZ2(),0.0,
            "Cll1122 or direct-MW contact entered Gmu matching",0.0,0.0);
    }
    // Both independent Cll contractions must still contribute to contact LL.
    const double expectedContact=c/(4.0*246.0*246.0);
    for(Cplx value:expected.contactAmplitudes(point()))
      close(value,expectedContact,"canonical LL contact has wrong normalization");
  }

  // All 2*2*4 mixed alias combinations retain the same physical strengths,
  // also in the presence of vertex operators and a coefficient-dependent MW.
  const SMEFTCalculator canonicalMixed(context,benchmark());
  for(int a=0;a<2;++a) for(int b=0;b<2;++b) for(int r=0;r<4;++r)
  {
    smeft::Input input=benchmark();
    input.setCll(one,one,two,two,0.0);
    input.setCll(one,two,two,one,0.0);
    input.setCee(one,one,two,two,0.0);
    input.setCll(ll[0][a][0],ll[0][a][1],ll[0][a][2],ll[0][a][3],.0004);
    input.setCll(ll[1][b][0],ll[1][b][1],ll[1][b][2],ll[1][b][3],.0001);
    input.setCee(rr[r][0],rr[r][1],rr[r][2],rr[r][3],.0002);
    sameContactPredictions(SMEFTCalculator(context,input),canonicalMixed,
                           "mixed contact alias population");
  }
}

void testExplicitRawContactImport(const ElectroweakSetup& context)
{
  // Historical Jackson raw entries carried independently added strengths.
  // Import their sums explicitly; aliases in the public API are assignments.
  const double ll1122=.0002, ll2211=-.0001, ll1221=.0003, ll2112=.00015;
  const double ee1122=.0001, ee1221=.0002, ee2211=-.00005, ee2112=.00015;
  const double a=ll1122+ll2211, b=ll1221+ll2112;
  const double rr=ee1122+ee1221+ee2211+ee2112;
  smeft::Input imported=zero(), alternative=zero();
  imported.setCll(one,one,two,two,a); imported.setCll(one,two,two,one,b);
  imported.setCee(one,one,two,two,rr);
  alternative.setCll(two,two,one,one,a); alternative.setCll(two,one,one,two,b);
  alternative.setCee(two,one,one,two,rr);
  const SMEFTCalculator expected(context,imported);
  sameContactPredictions(SMEFTCalculator(context,alternative),expected,
                         "explicit summed raw-contact import");
  const double scale=1.0/(4.0*246.0*246.0);
  const double signs[4]={1.0,-1.0,-1.0,1.0};
  const auto contacts=expected.contactAmplitudes(point());
  for(unsigned k=0;k<4;++k)
    close(contacts[k],scale*(a+b+signs[k]*rr),"raw-import summed contact strength");
  if(context.inputScheme()==EWInputScheme::GmuAlphaMZ)
    checkFiniteEquation(expected,smeft::WMassSolverOptions());

  smeft::Input assigned=zero();
  assigned.setCee(one,one,two,two,ee1122);
  assigned.setCee(one,two,two,one,ee1221);
  assigned.setCee(two,two,one,one,ee2211);
  assigned.setCee(two,one,one,two,ee2112);
  close(assigned.Cee(one,one,two,two),ee2112,
        "sequential alias setters did not use the latest value",0.0,0.0);
  require(assigned.Cee(one,one,two,two)!=rr,
          "raw import was incorrectly made implicit in alias setters");
}

SMEFTCalculator detachedPrediction()
{
  const SMval direct = makeSM();
  const SMvalGmu local(direct);
  const ElectroweakSetup context(EWInputScheme::GmuAlphaMZ,local);
  return SMEFTCalculator(context,benchmark());
}

void testOwnership()
{
  SMval direct = makeSM();
  smeft::Input coefficients = benchmark();
  const ElectroweakSetup context(EWInputScheme::AlphaMWMZ,direct);
  const SMEFTCalculator original(context,coefficients);
  SMEFTCalculator copied(original);
  SMEFTCalculator moved(std::move(copied));
  direct.set(MW,80.30); coefficients.clearWilsonCoefficients();
  samePrediction(original,moved,"copy/move or caller mutation changed snapshot");
  close(original.wilsonInput().CphiWB(),0.0002,"Wilson snapshot aliases caller");
  require(std::fabs(original.resolvedInput().get(InputPar::WMassComplex)
                    -direct.get(InputPar::WMassComplex))>0.01,
          "SM snapshot aliases caller");
  SMEFTCalculator assigned = detachedPrediction();
  assigned = original;
  samePrediction(original,assigned,"copy assignment changed owned snapshot");
  assigned = detachedPrediction();
  const SMEFTCalculator detached = detachedPrediction();
  samePrediction(assigned,detached,"temporary context or move assignment dangles");
  const double before = detached.differentialCrossSection(point()).total();
  for(unsigned i=0; i<4; ++i)
    close(detached.differentialCrossSection(point()).total(),before,
          "repeated call mutated native prediction",0.0,0.0);
}

void testValidation(const ElectroweakSetup& alpha, const ElectroweakSetup& gmu)
{
  const SMEFTCalculator prediction(alpha,zero());
  const double nan = std::numeric_limits<double>::quiet_NaN();
  const double inf = std::numeric_limits<double>::infinity();
  for(double invalid : {nan,inf,-1.0,0.0})
  {
    ScatteringPoint where = point(); where.s = invalid;
    throws([&]{ prediction.differentialCrossSection(where); },"invalid s accepted");
    where = point(); where.unitConversion = invalid;
    throws([&]{ prediction.amplitudes(where); },"invalid conversion accepted");
  }
  for(double invalid : {nan,inf,-1.0,1.0,1.01})
  {
    ScatteringPoint where = point(); where.cosTheta = invalid;
    throws([&]{ prediction.amplitudes(where); },"invalid/singular angle accepted");
  }
  throws([&]{ prediction.amplitudes(point(90.0,Fermion::electron)); },"Bhabha accepted");
  throws([&]{ prediction.amplitudes(point(90.0,Fermion::nuMuon)); },"neutrino process accepted");
  const double mz=prediction.resolvedInput().get(InputPar::ZMassComplex);
  throws([&]{ prediction.amplitudes(point(mz)); },"unsupported exact-pole background accepted");
  smeft::Input overflow = zero(1.0e-100);
  overflow.setCee(one,one,two,two,1.0e200);
  const SMEFTCalculator contactOverflow(alpha,overflow);
  throws([&]{ contactOverflow.contactAmplitudes(point()); },"nonfinite standalone contact accepted");
  throws([&]{ auto p=point(); p.initial=Fermion::muon; prediction.amplitudes(p); },
         "non-electron initial state accepted");
  throws([&]{ prediction.amplitude(point(),CurrentChannel(Current::scalar,Current::vector)); },
         "unsupported current accepted");
  throws([&]{ prediction.pseudoObservable(PseudoObservable::SW,Fermion::nuMuon); },
         "neutrino SW accepted");
  throws([&]{ prediction.pseudoObservable(PseudoObservable::DeltaR,Fermion::electron); },
         "undefined finite DeltaR accepted");
  // The private lookup must reject unsupported identifiers before any legacy
  // array is indexed, both through the calculator and standalone shifts.
  // Include negative/out-of-range indices and unsupported in-range species.
  const smeft::Input zeroCoefficients=zero();
  for(int invalid : {-1,0,6,10,17,20,999})
  {
    const Fermion fermion=static_cast<Fermion>(invalid);
    for(PseudoObservable observable : {PseudoObservable::FA,PseudoObservable::FV,
                                       PseudoObservable::SW})
      throwsInvalidArgument([&]{ prediction.pseudoObservable(observable,fermion); },
                            "invalid observable fermion accepted");
    throwsInvalidArgument([&]{ smeft::directZAxialShift(invalid,alpha.input(),zeroCoefficients); },
                          "invalid axial-shift fermion accepted");
    throwsInvalidArgument([&]{ smeft::directZVectorShift(invalid,alpha.input(),zeroCoefficients); },
                          "invalid vector-shift fermion accepted");
  }
  throws([&]{ SMEFTCalculator p(alpha,smeft::Input()); },"missing Lambda accepted");

  for(double invalid : {nan,inf,-1.0,0.0})
  {
    inval bad(alpha.input()); bad.set(InputPar::alpha,invalid);
    const ElectroweakSetup context(EWInputScheme::AlphaMWMZ,bad);
    throws([&]{ SMEFTCalculator p(context,zero()); },"invalid alpha accepted");
  }
  inval badMass(alpha.input()); badMass.set(InputPar::WMassComplex,100.0);
  const ElectroweakSetup badContext(EWInputScheme::AlphaMWMZ,badMass);
  throws([&]{ SMEFTCalculator p(badContext,zero()); },"MW >= MZ accepted");
  inval missing(alpha.input()); missing.set(InputPar::charmMass,nan);
  const ElectroweakSetup missingContext(EWInputScheme::AlphaMWMZ,missing);
  throws([&]{ SMEFTCalculator p(missingContext,zero()); },"missing high-SM input accepted");
  smeft::Input badFA = zero(); badFA.setCphil1(one,one,-10.0);
  throws([&]{ SMEFTCalculator p(alpha,badFA);
              p.pseudoObservable(PseudoObservable::FA,Fermion::electron); },
         "nonphysical negative FA accepted");

  smeft::WMassSolverOptions options;
  options.maxIterations=1;
  throws([&]{ SMEFTCalculator p(gmu,benchmark(),options); },"solver iteration limit ignored");
  options.maxIterations=0;
  throws([&]{ SMEFTCalculator p(gmu,zero(),options); },"zero solver iterations accepted");
  options = smeft::WMassSolverOptions(); options.massSquaredToleranceGeV2=nan;
  throws([&]{ SMEFTCalculator p(gmu,zero(),options); },"nonfinite solver tolerance accepted");
  options.massSquaredToleranceGeV2=0.0;
  throws([&]{ SMEFTCalculator p(gmu,zero(),options); },"zero solver tolerance accepted");
  smeft::Input huge = zero(); huge.setCphiWB(1000.0);
  throws([&]{ SMEFTCalculator p(gmu,huge); },"nonphysical mass solution accepted");
}

} // namespace

int main()
{
  try
  {
    const SMval direct = makeSM();
    const SMvalGmu derived(direct);
    const ElectroweakSetup alpha(EWInputScheme::AlphaMWMZ,direct);
    const ElectroweakSetup gmu(EWInputScheme::GmuAlphaMZ,derived);
    testZeroAndSquare(alpha); testZeroAndSquare(gmu);
    testFVAndLambda(alpha); testFVAndLambda(gmu);
    testFiniteGmu(gmu);
    testContactOnce(alpha); testContactOnce(gmu);
    testEquivalentContactLabels(alpha); testEquivalentContactLabels(gmu);
    testExplicitRawContactImport(alpha); testExplicitRawContactImport(gmu);
    testOwnership(); testValidation(alpha,gmu);
    std::cout << "Native SMEFT structural tests passed\n";
    return 0;
  }
  catch(const std::exception& error)
  {
    std::cerr << "Native SMEFT structural test failed: " << error.what() << '\n';
    return 1;
  }
}

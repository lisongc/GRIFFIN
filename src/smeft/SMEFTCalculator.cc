/*
  File: SMEFTCalculator.cc
  Purpose: Coordinate finite SMEFT and zero-coefficient SM evaluations and
           assemble scattering amplitudes/cross sections. Shared FA/SW/FV
           calculations are implemented in EWPO_smeft.cc.
  Authors: Jackson Wallace and Lisong Chen
  Last updated: 2026-11-16
*/

#include "griffin/smeft/SMEFTCalculator.h"

#include <cmath>
#include <stdexcept>
#include <string>

#include "griffin/sm/xscaas.h"
#include "griffin/sm/tools.h"

namespace griffin {
namespace prediction {
namespace {

// Fixed order of the public amplitude array: VV, AV, VA, AA.
const Current initialCurrents[4] = {
  Current::vector,Current::axial,Current::vector,Current::axial
};
const Current finalCurrents[4] = {
  Current::vector,Current::vector,Current::axial,Current::axial
};

smeft::Input zeroLike(const smeft::Input& coefficients)
{
  smeft::Input zero = coefficients;
  zero.clearWilsonCoefficients();
  return zero;
}

// Inspect a copy's protected storage before legacy get() can terminate the
// process for missing data. This does not modify or retain the caller's input.
class CheckedInput : public inval {
public:
  explicit CheckedInput(const inval& input) : inval(input) {}
  double require(InputPar key) const
  {
    const unsigned i = static_cast<unsigned>(key);
    if(i >= data.size() || !std::isfinite(data[i]))
      throw std::invalid_argument("SMEFT: missing/nonfinite SM input "
                                  +std::to_string(i));
    return data[i];
  }
};

void validateHighSMInput(const inval& input)
{
  const CheckedInput checked(input);
  const InputPar positive[] = {
    InputPar::alpha, InputPar::WMassComplex, InputPar::ZMassComplex,
    InputPar::HiggsMass, InputPar::topMass, InputPar::bottomMass,
    InputPar::ZWidthComplex
  };
  for(InputPar key : positive)
    if(checked.require(key) <= 0.0)
      throw std::domain_error("SMEFT: required positive SM input");
  const InputPar nonnegative[] = {
    InputPar::alphaS, InputPar::WWidthComplex, InputPar::electronMass,
    InputPar::muonMass, InputPar::tauMass, InputPar::downMass,
    InputPar::upMass, InputPar::strangeMass, InputPar::charmMass
  };
  for(InputPar key : nonnegative)
    if(checked.require(key) < 0.0)
      throw std::domain_error("SMEFT: negative SM mass/width/coupling");
  checked.require(InputPar::DeltaAlpha);
  if(checked.require(InputPar::WMassComplex)
      >= checked.require(InputPar::ZMassComplex))
    throw std::domain_error("SMEFT requires 0 < MW < MZ");
}

void validatePoint(const ScatteringPoint& point)
{
  if(!std::isfinite(point.s) || !std::isfinite(point.cosTheta)
      || !std::isfinite(point.unitConversion))
    throw std::invalid_argument("SMEFT cross-section kinematics must be finite");
  // Endpoint angles contain collinear singularities in the inherited virtual
  // background. This interface is not a real-radiation/inclusive observable.
  if(std::fabs(point.cosTheta) >= 1.0)
    throw std::domain_error("SMEFT cross sections require -1 < cos(theta) < 1");
  // Validate before constructing amplitudes. The conversion utility also
  // checks its own inputs, but must not be used as a dummy calculation here.
  if(point.initial != Fermion::electron)
    throw std::invalid_argument("SMEFT scattering requires an electron initial state");
  switch(point.final)
  {
    case Fermion::muon: case Fermion::tau:
    case Fermion::d: case Fermion::u: case Fermion::s:
    case Fermion::c: case Fermion::b:
      break;
    case Fermion::electron:
      throw std::domain_error("Bhabha scattering requires dedicated s/t-channel structures");
    case Fermion::nuElectron: case Fermion::nuMuon: case Fermion::nuTau:
      throw std::domain_error("neutrino final states require a dedicated process implementation");
    default:
      throw std::invalid_argument("unsupported final state for SMEFT scattering");
  }
  if(point.s <= 0.0)
    throw std::domain_error("cross-section invariant mass squared must be positive");
  if(point.unitConversion <= 0.0)
    throw std::domain_error("cross-section unit conversion must be positive");
}

class FixedObservable : public psobs {
  Cplx value_;
public:
  FixedObservable(Cplx value, const inval& input) : psobs(input), value_(value) {}
  Cplx result() const override { return value_; }
};

} // namespace

void SMEFTCalculator::requireFinite(Cplx value, const char* quantity)
{
  if(!std::isfinite(real(value)) || !std::isfinite(imag(value)))
    throw std::domain_error(std::string("nonfinite SMEFT ")+quantity);
}

SMEFTCalculator::Amplitudes SMEFTCalculator::amplitudesAt(
    const ScatteringPoint& point, const inval& input,
    const smeft::Input& coefficients)
{
  validatePoint(point);
  const double mz = input.get(InputPar::ZMassComplex);
  if(point.s == mz*mz)
    throw std::domain_error("exact-real-Z-pole evaluation is unsupported; "
                            "near-pole precision requires dedicated validation. "
                            "See docs/KNOWN_ISSUES.md");
  const PoleInputs initial = poleInputs(point.initial,input,coefficients);
  const PoleInputs final = poleInputs(point.final,input,coefficients);
  const FixedObservable fai(initial.fa,input), fao(final.fa,input);
  const FixedObservable swi(initial.sw,input), swo(final.sw,input);
  SMEFTCalculator::Amplitudes result;
  for(std::size_t k = 0; k < result.size(); ++k)
  {
    const Current ci = initialCurrents[k];
    const Current cf = finalCurrents[k];
    // Include the selected mixed electroweak-QCD residue and off-resonance terms at
    // this resolved mass, for both the full prediction and its SM result.
    // The matrix element borrows the four local observables during evaluation.
    const mat_SMaas poleAndBackground(point.initial,point.final,ci,cf,
        fai,fao,swi,swo,point.s,point.cosTheta,input);
    result[k] = poleAndBackground.result()+smeft::fourFermionContactAmplitude(
        point.initial,point.final,ci,cf,coefficients);
    requireFinite(result[k],"amplitude");
  }
  return result;
}

SMEFTCalculator::SMEFTCalculator(
    const ElectroweakSetup& setup, const smeft::Input& coefficients,
    const smeft::WMassSolverOptions& options)
  : scheme_(setup.inputScheme()), coefficients_(coefficients),
    zeroCoefficients_(zeroLike(coefficients)),
    // Resolve the SM input values separately with zero coefficients; they are
    // not the shifted inputs with only direct vertex/contact terms switched off.
    smMassSolution_(setup.input(),zeroCoefficients_,scheme_,options),
    point_(setup.input(),coefficients_,scheme_,options)
{
  validateHighSMInput(smInputValues());
  validateHighSMInput(resolvedInput());
}

FinitePrediction<Cplx> SMEFTCalculator::pseudoObservable(
    PseudoObservable observable, Fermion fermion) const
{
  return FinitePrediction<Cplx>(
      observableAt(observable,fermion,smInputValues(),zeroCoefficients_),
      observableAt(observable,fermion,resolvedInput(),coefficients_));
}

SMEFTCalculator::Amplitudes SMEFTCalculator::amplitudes(
    const ScatteringPoint& point) const
{
  return amplitudesAt(point,resolvedInput(),coefficients_);
}

SMEFTCalculator::Amplitudes SMEFTCalculator::smAmplitudes(
    const ScatteringPoint& point) const
{
  return amplitudesAt(point,smInputValues(),zeroCoefficients_);
}

SMEFTCalculator::Amplitudes SMEFTCalculator::contactAmplitudes(
    const ScatteringPoint& point) const
{
  validatePoint(point);
  Amplitudes result;
  for(std::size_t k = 0; k < result.size(); ++k)
  {
    result[k] = smeft::fourFermionContactAmplitude(point.initial,point.final,
        initialCurrents[k],finalCurrents[k],coefficients_);
    requireFinite(result[k],"contact amplitude");
  }
  return result;
}

FinitePrediction<Cplx> SMEFTCalculator::amplitude(
    const ScatteringPoint& point, const CurrentChannel& channel) const
{
  if((channel.initial != Current::vector && channel.initial != Current::axial)
      || (channel.final != Current::vector && channel.final != Current::axial))
    throw std::domain_error("SMEFT supports only vector/axial currents");
  const unsigned k = (channel.initial == Current::axial ? 1 : 0)
                     +(channel.final == Current::axial ? 2 : 0);
  return FinitePrediction<Cplx>(smAmplitudes(point)[k],amplitudes(point)[k]);
}

FinitePrediction<double> SMEFTCalculator::differentialCrossSection(
    const ScatteringPoint& point) const
{
  validatePoint(point);
  // Square each complete amplitude. Their finite difference includes the
  // nonlinear terms generated by this prescription, not just SM-D6 interference.
  const double smResult = differentialCrossSectionFromAmplitudes(
      smAmplitudes(point),point.initial,point.final,
      point.s,point.cosTheta,point.unitConversion);
  const double total = differentialCrossSectionFromAmplitudes(
      amplitudes(point),point.initial,point.final,
      point.s,point.cosTheta,point.unitConversion);
  if(!std::isfinite(smResult) || !std::isfinite(total))
    throw std::domain_error("nonfinite SMEFT cross section");
  return FinitePrediction<double>(smResult,total);
}

} // namespace prediction
} // namespace griffin

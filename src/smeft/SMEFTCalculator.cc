#include "griffin/SMEFTCalculator.h"

#include <cmath>
#include <stdexcept>
#include <string>

#include "griffin/sm/EWPOZ2.h"
#include "griffin/sm/ff.h"
#include "griffin/sm/ff0.h"
#include "smeft/detail/SMEFTFlavor.h"
#include "griffin/smeft/ff0_smeft.h"
#include "griffin/smeft/xsc_smeft.h"
#include "griffin/sm/xscaas.h"

namespace griffin {
namespace prediction {
namespace {

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

void requireFinite(Cplx value, const char* quantity)
{
  if(!std::isfinite(real(value)) || !std::isfinite(imag(value)))
    throw std::domain_error(std::string("nonfinite SMEFT ")+quantity);
}

struct PoleInputs { Cplx fa, sw; };

PoleInputs poleInputs(Fermion fermion, const inval& input,
                      const smeft::Input& coefficients)
{
  const int type = index(fermion);
  smeft::detail::fermionInfo(type); // validate before indexing legacy arrays
  const double a0 = az0(type,input);
  const double v0 = vz0(type,input);
  // The input already contains the resolved MW, including its finite Gmu
  // shift.  These helpers add only the explicit vertex changes at that point.
  const double da = smeft::directZAxialShift(type,input,coefficients);
  const double dv = smeft::directZVectorShift(type,input,coefficients);
  PoleInputs result;
  // These explicit FA/SW additions are linear in Wilson coefficients at a
  // fixed mass point. In particular FA has no explicit da*da term. The SM
  // functions are NNLO+ and use the resolved mass. Downstream pole/FV
  // construction and full amplitude squaring retain their generated nonlinear
  // terms; this is not a prescription for all higher powers of the coefficients.
  result.fa = FA_SMNNLO(type,input).result()+2.0*a0*da;
  result.sw = 0.0; // unused for neutrinos; their public SW request is rejected
  if(!isNeutrino(fermion))
    result.sw = SW_SMNNLO(type,input).result()
      -(dv/a0-v0*da/(a0*a0))/(4.0*std::fabs(Qf[type]));
  requireFinite(result.fa,"FA");
  requireFinite(result.sw,"SW");
  if(real(result.fa) <= 0.0)
    throw std::domain_error("SMEFT FA must be positive for the pole construction");
  return result;
}

Cplx observableAt(PseudoObservable observable, Fermion fermion,
                   const inval& input, const smeft::Input& coefficients)
{
  if(observable == PseudoObservable::WMass)
    return input.get(InputPar::WMassComplex);
  if(observable == PseudoObservable::DeltaR)
    throw std::domain_error("finite DeltaR is not defined by SMEFTCalculator; "
                            "use massSolution() for the Gmu equation diagnostics");
  if(observable == PseudoObservable::SW && isNeutrino(fermion))
    throw std::domain_error("SW is undefined for a neutral fermion");
  const PoleInputs pole = poleInputs(fermion,input,coefficients);
  if(observable == PseudoObservable::FA) return pole.fa;
  if(observable == PseudoObservable::SW) return pole.sw;
  if(observable != PseudoObservable::FV)
    throw std::invalid_argument("unsupported SMEFT pseudo-observable");

  // Use the FV_SMNNLO construction with the same shifted FA and SW
  // as the pole amplitude. Retain the SM imaginary-ratio term at resolved MW;
  // otherwise even the zero-coefficient limit would differ from GRIFFIN SM.
  const int type = index(fermion);
  const double a0 = az0(type,input), v0 = vz0(type,input);
  const double imaginaryRatio =
      (a0*(iz1f(type,VEC,input)+iz1b(type,VEC,input))
       -v0*(iz1f(type,AXV,input)+iz1b(type,AXV,input)))/(a0*a0);
  const double q = 1.0-4.0*std::fabs(Qf[type])*real(pole.sw);
  const Cplx fv = pole.fa*(q*q+imaginaryRatio*imaginaryRatio);
  requireFinite(fv,"FV");
  return fv;
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
  const FermionPairCrossSection validate(
      point.initial,point.final,point.s,point.cosTheta,point.unitConversion);
  (void)validate;
}

class FixedObservable : public psobs {
  Cplx value_;
public:
  FixedObservable(Cplx value, const inval& input) : psobs(input), value_(value) {}
  Cplx result() const override { return value_; }
};

SMEFTCalculator::Amplitudes amplitudesAt(
    const ScatteringPoint& point, const inval& input,
    const smeft::Input& coefficients)
{
  validatePoint(point);
  const double mz = input.get(InputPar::ZMassComplex);
  if(point.s == mz*mz)
    throw std::domain_error("exact-real-Z-pole evaluation is unsupported; nearby "
                            "mat_SMaas scattering also has unresolved limitations. "
                            "See docs/KNOWN_ISSUES.md");
  const PoleInputs initial = poleInputs(point.initial,input,coefficients);
  const PoleInputs final = poleInputs(point.final,input,coefficients);
  const FixedObservable fai(initial.fa,input), fao(final.fa,input);
  const FixedObservable swi(initial.sw,input), swo(final.sw,input);
  SMEFTCalculator::Amplitudes result;
  for(std::size_t k = 0; k < result.size(); ++k)
  {
    const Current ci = FermionPairCrossSection::initialCurrent(k);
    const Current cf = FermionPairCrossSection::finalCurrent(k);
    // Include the selected mixed electroweak-QCD residue and off-resonance terms at
    // this resolved mass, for both the full prediction and its SM reference.
    // The borrowed-observable constructor avoids legacy raw-pointer ownership.
    const mat_SMaas poleAndBackground(point.initial,point.final,ci,cf,
        fai,fao,swi,swo,point.s,point.cosTheta,input);
    result[k] = poleAndBackground.result()+smeft::coeffS4fD6Tree(
        point.initial,point.final,ci,cf,coefficients);
    requireFinite(result[k],"amplitude");
  }
  return result;
}

} // namespace

SMEFTCalculator::SMEFTCalculator(
    const ElectroweakSetup& setup, const smeft::Input& coefficients,
    const smeft::WMassSolverOptions& options)
  : scheme_(setup.inputScheme()), coefficients_(coefficients),
    zeroCoefficients_(zeroLike(coefficients)),
    // Resolve the zero-coefficient point separately; it is not the shifted
    // point with only its explicit vertex/contact terms switched off.
    referencePoint_(setup.input(),zeroCoefficients_,scheme_,options),
    point_(setup.input(),coefficients_,scheme_,options)
{
  validateHighSMInput(referenceInput());
  validateHighSMInput(resolvedInput());
}

FinitePrediction<Cplx> SMEFTCalculator::pseudoObservable(
    PseudoObservable observable, Fermion fermion) const
{
  return FinitePrediction<Cplx>(
      observableAt(observable,fermion,referenceInput(),zeroCoefficients_),
      observableAt(observable,fermion,resolvedInput(),coefficients_));
}

SMEFTCalculator::Amplitudes SMEFTCalculator::amplitudes(
    const ScatteringPoint& point) const
{
  return amplitudesAt(point,resolvedInput(),coefficients_);
}

SMEFTCalculator::Amplitudes SMEFTCalculator::referenceAmplitudes(
    const ScatteringPoint& point) const
{
  return amplitudesAt(point,referenceInput(),zeroCoefficients_);
}

SMEFTCalculator::Amplitudes SMEFTCalculator::contactAmplitudes(
    const ScatteringPoint& point) const
{
  validatePoint(point);
  Amplitudes result;
  for(std::size_t k = 0; k < result.size(); ++k)
  {
    result[k] = smeft::coeffS4fD6Tree(point.initial,point.final,
        FermionPairCrossSection::initialCurrent(k),
        FermionPairCrossSection::finalCurrent(k),coefficients_);
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
  return FinitePrediction<Cplx>(referenceAmplitudes(point)[k],amplitudes(point)[k]);
}

FinitePrediction<double> SMEFTCalculator::differentialCrossSection(
    const ScatteringPoint& point) const
{
  validatePoint(point);
  const FermionPairCrossSection crossSection(
      point.initial,point.final,point.s,point.cosTheta,point.unitConversion);
  // Square each complete amplitude. Their finite difference includes the
  // nonlinear terms generated by this prescription, not just SM-D6 interference.
  const double reference = crossSection.squaredPrediction(referenceAmplitudes(point));
  const double total = crossSection.squaredPrediction(amplitudes(point));
  if(!std::isfinite(reference) || !std::isfinite(total))
    throw std::domain_error("nonfinite SMEFT cross section");
  return FinitePrediction<double>(reference,total);
}

} // namespace prediction
} // namespace griffin

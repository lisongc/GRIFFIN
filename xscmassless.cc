#include "xscmassless.h"

#include <stdexcept>

namespace griffin {
namespace {

bool isSupportedChargedFinalState(Fermion fermion)
{
  switch(fermion)
  {
    case Fermion::d:
    case Fermion::u:
    case Fermion::s:
    case Fermion::c:
    case Fermion::b:
    case Fermion::muon:
    case Fermion::tau:
      return true;
    default:
      return false;
  }
}

double finalStateColorMultiplicity(Fermion fermion)
{
  switch(fermion)
  {
    case Fermion::d:
    case Fermion::u:
    case Fermion::s:
    case Fermion::c:
    case Fermion::b:
      return 3.0;
    default:
      return 1.0;
  }
}

double angularSquare(
    const MasslessFermionCrossSectionBuilder::Amplitudes& amplitude,
    double cosTheta)
{
  const double even = 1.0+cosTheta*cosTheta;
  const double odd = 4.0*cosTheta;
  const Cplx squared = even*(
      amplitude[0]*conj(amplitude[0])
      + amplitude[1]*conj(amplitude[1])
      + amplitude[2]*conj(amplitude[2])
      + amplitude[3]*conj(amplitude[3]))
      + odd*(amplitude[0]*conj(amplitude[3])
             + amplitude[2]*conj(amplitude[1]));
  return real(squared);
}

double angularInterference(
    const MasslessFermionCrossSectionBuilder::Amplitudes& reference,
    const MasslessFermionCrossSectionBuilder::Amplitudes& correction,
    double cosTheta)
{
  const double even = 1.0+cosTheta*cosTheta;
  const double odd = 4.0*cosTheta;
  const Cplx interference = even*(
      correction[0]*conj(reference[0])
      + reference[0]*conj(correction[0])
      + correction[1]*conj(reference[1])
      + reference[1]*conj(correction[1])
      + correction[2]*conj(reference[2])
      + reference[2]*conj(correction[2])
      + correction[3]*conj(reference[3])
      + reference[3]*conj(correction[3]))
      + odd*(correction[0]*conj(reference[3])
             + reference[0]*conj(correction[3])
             + correction[2]*conj(reference[1])
             + reference[2]*conj(correction[1]));
  return real(interference);
}

} // namespace

MasslessFermionCrossSectionBuilder::MasslessFermionCrossSectionBuilder(
    Fermion initial, Fermion final, double invariantMassSquared,
    double cosine, double unitConversion)
  : initial_(initial), final_(final), s_(invariantMassSquared),
    cosTheta_(cosine), unitConversion_(unitConversion)
{
  if(initial_ != Fermion::electron)
    throw std::invalid_argument(
        "massless annihilation cross section currently requires an electron initial state");
  if(final_ == Fermion::electron)
    throw std::domain_error(
        "Bhabha scattering requires dedicated s/t-channel and scalar/pseudoscalar structures");
  if(isNeutrino(final_))
    throw std::domain_error(
        "neutrino final states require a dedicated process implementation");
  if(!isSupportedChargedFinalState(final_))
    throw std::invalid_argument(
        "unsupported final state for the massless annihilation cross section");
  if(s_ <= 0.0 || cosTheta_ < -1.0 || cosTheta_ > 1.0)
    throw std::domain_error(
        "cross-section kinematics are outside the physical domain");
  if(unitConversion_ <= 0.0)
    throw std::domain_error("cross-section unit conversion must be positive");
}

Current MasslessFermionCrossSectionBuilder::initialCurrent(
    std::size_t channel)
{
  static const Current currents[4] = {
    Current::vector,Current::axial,Current::vector,Current::axial
  };
  if(channel >= channelCount())
    throw std::out_of_range("massless current-channel index is out of range");
  return currents[channel];
}

Current MasslessFermionCrossSectionBuilder::finalCurrent(
    std::size_t channel)
{
  static const Current currents[4] = {
    Current::vector,Current::vector,Current::axial,Current::axial
  };
  if(channel >= channelCount())
    throw std::out_of_range("massless current-channel index is out of range");
  return currents[channel];
}

double MasslessFermionCrossSectionBuilder::squaredPrediction(
    const Amplitudes& amplitudes) const
{
  const double prefactor = finalStateColorMultiplicity(final_)*s_
      /(32.0*Pi)*unitConversion_;
  return prefactor*angularSquare(amplitudes,cosTheta_);
}

double MasslessFermionCrossSectionBuilder::linearInterference(
    const Amplitudes& reference, const Amplitudes& correction) const
{
  const double prefactor = finalStateColorMultiplicity(final_)*s_
      /(32.0*Pi)*unitConversion_;
  return prefactor*angularInterference(reference,correction,cosTheta_);
}

} // namespace griffin

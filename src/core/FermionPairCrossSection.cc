#include "griffin/core/FermionPairCrossSection.h"

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
    const FermionPairCrossSection::Amplitudes& amplitude,
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

} // namespace

FermionPairCrossSection::FermionPairCrossSection(
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

Current FermionPairCrossSection::initialCurrent(
    std::size_t channel)
{
  static const Current currents[4] = {
    Current::vector,Current::axial,Current::vector,Current::axial
  };
  if(channel >= channelCount())
    throw std::out_of_range("massless current-channel index is out of range");
  return currents[channel];
}

Current FermionPairCrossSection::finalCurrent(
    std::size_t channel)
{
  static const Current currents[4] = {
    Current::vector,Current::vector,Current::axial,Current::axial
  };
  if(channel >= channelCount())
    throw std::out_of_range("massless current-channel index is out of range");
  return currents[channel];
}

double FermionPairCrossSection::squaredPrediction(
    const Amplitudes& amplitudes) const
{
  const double prefactor = finalStateColorMultiplicity(final_)*s_
      /(32.0*Pi)*unitConversion_;
  return prefactor*angularSquare(amplitudes,cosTheta_);
}

} // namespace griffin

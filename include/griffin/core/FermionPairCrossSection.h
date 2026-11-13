/*
  FermionPairCrossSection.h: differential cross section from fermion-pair amplitudes

  Combine the supplied VV, AV, VA and AA amplitudes into unpolarized
  d sigma / d cos(theta), treating external fermions as massless.
  Currently supports e+e- annihilation to muons, taus and u,d,s,c,b;
  Bhabha scattering and neutrino final states require separate implementations.
  This helper does not compute SM/SMEFT amplitudes or choose their accuracy.
*/

#ifndef GRIFFIN_FERMION_PAIR_CROSS_SECTION_H
#define GRIFFIN_FERMION_PAIR_CROSS_SECTION_H

#include <array>
#include <cstddef>

#include "griffin/core/classes.h"

namespace griffin {

class FermionPairCrossSection {
public:
  typedef std::array<Cplx,4> Amplitudes;

private:
  Fermion initial_;
  Fermion final_;
  double s_;
  double cosTheta_;
  double unitConversion_;

public:
  FermionPairCrossSection(
      Fermion initial, Fermion final, double invariantMassSquared,
      double cosine, double unitConversion = 1.0);

  static std::size_t channelCount() { return 4; }
  static Current initialCurrent(std::size_t channel);
  static Current finalCurrent(std::size_t channel);

  Fermion initial() const { return initial_; }
  Fermion final() const { return final_; }
  double invariantMassSquared() const { return s_; }
  double cosTheta() const { return cosTheta_; }
  double unitConversion() const { return unitConversion_; }

  double squaredPrediction(const Amplitudes& amplitudes) const;
};

} // namespace griffin

#endif // GRIFFIN_FERMION_PAIR_CROSS_SECTION_H

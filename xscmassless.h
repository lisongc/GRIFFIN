/*
  xscmassless.h: shared massless-fermion differential cross-section builder

  This builder owns the process-domain contract and the common vector/axial
  angular algebra.  Model-specific code supplies the four current amplitudes;
  it must not duplicate the observable-level square or interference formula.
*/

#ifndef GRIFFIN_XSCMASSLESS_H
#define GRIFFIN_XSCMASSLESS_H

#include <array>
#include <cstddef>

#include "classes.h"

namespace griffin {

class MasslessFermionCrossSectionBuilder {
public:
  typedef std::array<Cplx,4> Amplitudes;

private:
  Fermion initial_;
  Fermion final_;
  double s_;
  double cosTheta_;
  double unitConversion_;

public:
  MasslessFermionCrossSectionBuilder(
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
  double linearInterference(const Amplitudes& reference,
                            const Amplitudes& correction) const;
};

} // namespace griffin

#endif // GRIFFIN_XSCMASSLESS_H

/*
  File: xsc_smeft.h
  Purpose: Define scattering kinematics/current requests for SMEFTCalculator
           and declare local dimension-six four-fermion contact amplitudes.
           The reduced contact amplitude is independent of s. SMEFTCalculator
           adds it once before taking the full amplitude square.
  Authors: Jackson Wallace and Lisong Chen
  Last updated: 2026-11-16
*/

#ifndef GRIFFIN_SMEFT_XSC_H
#define GRIFFIN_SMEFT_XSC_H

#include "griffin/core/classes.h"
#include "griffin/smeft/SMEFTInput.h"

namespace griffin {
namespace prediction {

struct ScatteringPoint {
  Fermion initial;       // fermion in the incoming fermion-antifermion pair
  Fermion final;         // fermion in the outgoing fermion-antifermion pair
  double s;              // center-of-mass energy squared, in GeV^2
  double cosTheta;       // cosine of the scattering angle, not the angle
  double unitConversion; // multiply GeV^-2 by this factor (0.38937966e9 for pb)

  ScatteringPoint(Fermion initialState, Fermion finalState,
                  double invariantMassSquared, double cosine,
                  double conversion = 1.0)
    : initial(initialState), final(finalState),
      s(invariantMassSquared), cosTheta(cosine),
      unitConversion(conversion) {}
};

struct CurrentChannel {
  Current initial;
  Current final;

  CurrentChannel(Current initialCurrent, Current finalCurrent)
    : initial(initialCurrent), final(finalCurrent) {}
};

} // namespace prediction

namespace smeft {

/*
  Massless vector/axial matching documented in docs/SMEFT_IMPLEMENTATION.md.
  Electron initial states only; Bhabha, scalar and pseudoscalar channels
  require separate process structures and are rejected. Neutrino contacts
  are available as coefficients, not as a complete scattering prediction.
*/
Cplx fourFermionContactAmplitude(int intype, int outtype, int inform, int outform,
                                const Input& eftInput);

inline Cplx fourFermionContactAmplitude(Fermion intype, Fermion outtype,
                                       Current inform, Current outform,
                                       const Input& eftInput)
{
  return fourFermionContactAmplitude(index(intype),index(outtype),index(inform),
                                    index(outform),eftInput);
}

} // namespace smeft
} // namespace griffin

#endif // GRIFFIN_SMEFT_XSC_H

/*
  Standalone local four-fermion contact amplitudes for SMEFTCalculator.

  Their reduced amplitude is independent of s, so its S coefficient is the
  contact amplitude itself and its S' coefficient vanishes. These functions
  neither solve for MW nor square the amplitude. The D6Tree suffix describes
  these tree-level dimension-six contact amplitudes, not an approximation
  order imposed on the complete SMEFTCalculator prediction.
*/

#ifndef GRIFFIN_SMEFT_XSC_H
#define GRIFFIN_SMEFT_XSC_H

#include "griffin/core/classes.h"
#include "griffin/smeft/SMEFTInput.h"

namespace griffin {
namespace smeft {

/*
  Massless vector/axial matching documented in docs/SMEFT_IMPLEMENTATION.md.
  Electron initial states only; Bhabha, scalar and pseudoscalar channels
  require separate process structures and are rejected. Neutrino contacts
  are available as coefficients, not as a complete scattering prediction.
*/
Cplx coeffS4fD6Tree(int intype, int outtype, int inform, int outform,
                    const Input& eftInput);
Cplx coeffSp4fD6Tree(int intype, int outtype, int inform, int outform,
                     const Input& eftInput);

inline Cplx coeffS4fD6Tree(Fermion intype, Fermion outtype, Current inform,
                           Current outform, const Input& eftInput)
{
  return coeffS4fD6Tree(index(intype),index(outtype),index(inform),
                        index(outform),eftInput);
}

inline Cplx coeffSp4fD6Tree(Fermion intype, Fermion outtype, Current inform,
                            Current outform, const Input& eftInput)
{
  return coeffSp4fD6Tree(index(intype),index(outtype),index(inform),
                         index(outform),eftInput);
}

} // namespace smeft
} // namespace griffin

#endif // GRIFFIN_SMEFT_XSC_H

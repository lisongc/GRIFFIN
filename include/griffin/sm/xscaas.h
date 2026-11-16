#ifndef GRIFFIN_SM_XSCAAS_H
#define GRIFFIN_SM_XSCAAS_H

/*
  File: xscaas.h
  Purpose: Declare the mixed electroweak-QCD extension of mat_SMNNLO.
  Upstream author: Ayres Freitas
  Merged integration: Jackson Wallace and Lisong Chen
  Last updated: 2026-11-16
*/

#include "griffin/sm/xscnnlo.h"

namespace griffin {

// SM matrix element including the selected mixed fermionic electroweak-QCD
// O(alpha_f alpha_s) extension in the residue and off-resonance remainder.
// This extends mat_SMNNLO; it does not supply complete NNLO off-resonance
// corrections. The added remainder supports different initial/final flavors only.
// Includes the upstream v1.1.1 corrected grids and near-pole Z self-energy
// interpolation. Literal evaluation at s=MZ^2 still needs a separate limit;
// SMEFTCalculator rejects that point. See docs/KNOWN_ISSUES.md for tested scope.
class mat_SMaas : public mat_SMNNLO {
public:
  using mat_SMNNLO::mat_SMNNLO;
  Cplx coeffR(void) const;     // complete residue, including the extension
  Cplx resoffZaas(void) const; // O(af as) off-resonance contribution
  Cplx resoffZ(void) const     // total off-resonance contribution
  {
    return(mat_SMNNLO::resoffZ() + resoffZaas());
  }
  Cplx result(void) const;
};

} // namespace

#endif // GRIFFIN_SM_XSCAAS_H

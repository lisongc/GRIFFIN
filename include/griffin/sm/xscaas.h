#ifndef GRIFFIN_SM_XSCAAS_H
#define GRIFFIN_SM_XSCAAS_H

/* xscaas.h: header file for xscaas.cc */

#include "griffin/sm/xscnnlo.h"

namespace griffin {

// SM matrix element including the selected mixed fermionic electroweak-QCD
// O(alpha_f alpha_s) extension in the residue and off-resonance remainder.
// This extends mat_SMNNLO; it does not supply complete NNLO off-resonance
// corrections. The added remainder supports different initial/final flavors only.
// Includes upstream v1.1.1 mixed-self-energy stabilization and v1.2 one-loop
// interpolation. The exact real Z pole remains unsupported. Near-pole
// precision and interpolation boundaries need dedicated validation; see
// docs/KNOWN_ISSUES.md and the background diagnostic.
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

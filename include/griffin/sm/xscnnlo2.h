#ifndef GRIFFIN_SM_XSCNNLO2_H
#define GRIFFIN_SM_XSCNNLO2_H

/* xscnnlo2.h: header file for xscnnlo2.cc */

#include "griffin/sm/xscaas.h"

namespace griffin {

// Selected fermionic two-loop off-peak additions to mat_SMaas, imported from
// upstream ca2d523. The upstream off-shell fermion-boson vertex/self-energy
// term and pure-photon vacuum-polarization term are currently disabled.
// This low-level class is separate from SMEFTCalculator's fixed mat_SMaas path.
// It requires positive charged-fermion masses. Its fermion-boson sector is
// incomplete for bottom, neutrino and identical-flavor channels, and the new
// box contribution is not validated near the Z pole. See docs/KNOWN_ISSUES.md.
class mat_SMNNLO2 : public mat_SMaas {
private:
  Cplx compoffZ2f(double sval) const;  // for computation of O(af^2) off-resonance contribution
  Cplx compoffZfb(double sval) const;  // for computation of O(af*fb) off-resonance contribution
  Cplx compboxfb(double sval) const;   // for computation of O(af*fb) box contribution
public:
  using mat_SMaas::mat_SMaas;
  Cplx resoffZ2f(void) const;  // O(af^2) off-resonance contribution
  Cplx resoffZfb(void) const;  // O(af*ab) off-resonance contribution
  Cplx resoffZ(void) const     // total
  {
    return(mat_SMaas::resoffZ()+resoffZ2f()+resoffZfb());
  }
  Cplx resvpol(void) const;	// added by result(), separately from resoffZ()
  Cplx result(void) const;
};

} // namespace

#endif // GRIFFIN_SM_XSCNNLO2_H

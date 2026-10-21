/* xscnnlo2.h: header file for xscnnlo2.cc */

#include "xscaas.h"

namespace griffin {

// NNLO matrix element predicted in the SM, including off-peak
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
  Cplx resvpol(void) const;	// contributions from photon vacuum polarization
  Cplx result(void) const;
};

} // namespace

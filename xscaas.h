/* xscaas.h: header file for xscaas.cc */

#include "xscnnlo.h"

namespace griffin {

// SM matrix element with NNLO O(af as)-corrections also off-peak
class mat_SMaas : public mat_SMNNLO {
public:
  using mat_SMNNLO::mat_SMNNLO;
  Cplx coeffR(void) const;     // correction to R coefficient
  Cplx resoffZaas(void) const; // O(af as) off-resonance contribution
  Cplx resoffZ(void) const     // total off-resonance contribution
  {
    return(mat_SMNNLO::resoffZ() + resoffZaas());
  }
  Cplx result(void) const;
};

} // namespace

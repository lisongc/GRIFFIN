/* xscSMEFT.h: header file for xscSMEFT.cc */

#ifndef __xscSMEFT__
#define __xscSMEFT__

#include "xscnnlo.h"

namespace griffin {

// matrix element SMEFT (interference term off-peak)
class mat_SMEFTLO : public mat_SMNNLO {
public:
  using mat_SMNNLO::mat_SMNNLO;

  // corrections to off-resonance contribution: 
  double resoff4f(void) const;  // four-fermion correction
  Cplx resoffZ(void) const     // total
  {
    return(mat_SMNNLO::resoffZ()+resoff4f());
  }
  Cplx result(void) const;
};

} // namespace

#endif // __xscSMEFT__

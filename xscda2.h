/* xscda2.h: header file for xscda2.cc */

#include "xscnnlo.h"

namespace griffin {

// matrix element predicted in the SM with Delal^2 terms in photon amplitude
class mat_SMda2 : public mat_SMNNLO {
public:
  using mat_SMNNLO::mat_SMNNLO;
  Cplx resoffZ2f(void) const;
  Cplx resoffZ(void) const
  {
    return(mat_SMNNLO::resoffZ()+resoffZ2f());
  }
  Cplx result(void) const;
};

} // namespace

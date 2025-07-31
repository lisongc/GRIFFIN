/* xscrunal.h: header file for xscrunal.cc */

#include "xscnnlo.h"

namespace griffin {

#define alQs 92
#define alQt 93

// matrix element predicted as in mat_SMNNLO, but with running alpha(s/t) 
// in the s/t-channel QED amplitudes
class mat_SMNNLOrunal : public mat_SMNNLO {
public:
  using mat_SMNNLO::mat_SMNNLO;

  // corrections to off-resonance contribution: 
  Cplx resoffZ0(void) const;   // tree-level
  Cplx resoffZ1f(void) const;  // 1-loop with closed fermion loops
  Cplx resoffZ(void) const     // total
  {
    return(resoffZ0()+resoffZ1f()+mat_SMNNLO::resoffZ1b());
  }
  Cplx result(void) const;
};

} // namespace

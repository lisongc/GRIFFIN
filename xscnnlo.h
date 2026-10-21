/* xscnnlo.h: header file for xscnnlo.cc */

#ifndef __xscnnlo__
#define __xscnnlo__

#include "classes.h"

namespace griffin {

// matel @ NLO is same as matel @ LO
#define mat_SMNLO matel

// matrix element predicted in the SM (NNLO on Z peak, NLO off-peak)
class mat_SMNNLO : public matel {
public:
  using matel::matel;
  Cplx coeffR(void) const;    // correction to R coefficient
  
  // corrections to S coefficient:
  Cplx coeffS1f(void) const;  // 1-loop with closed fermion loops
  Cplx coeffS1b(void) const;  // 1-loop without closed fermion loops 
  Cplx coeffS(void) const     // total
  {
    return(matel::coeffS()+coeffS1f()+coeffS1b());
  }

  // corrections to off-resonance contribution: 
  Cplx resoffZ1f(void) const;  // 1-loop with closed fermion loops
  Cplx resoffZ1b(void) const;  // 1-loop without closed fermion loops
  Cplx resoffZ(void) const     // total
  {
    return(matel::resoffZ()+resoffZ1f()+resoffZ1b());
  }
  Cplx result(void) const;
};

} // namespace

#endif // __xscnnlo__

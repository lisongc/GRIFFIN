/* xscenlo.h: header file for xscenlo.cc */

#include "classes.h"
#include "ff0.h"

namespace griffin {

// matrix element at EXACTLY NLO (no products of 1-loop terms), for testing
class mat_SMeNLO : public matel {
public:
  using matel::matel;
  Cplx coeffR(void) const;    // correction to R coefficient
  
  // corrections to off-resonance contribution: 
  Cplx resoffZ1f(void) const;  // 1-loop with closed fermion loops
  Cplx resoffZ1b(void) const;  // 1-loop without closed fermion loops
  Cplx resoffZ(void) const     // total
  {
    return(g0(it,iff,*ival)*g0(ot,off,*ival)/s+resoffZ1f()+resoffZ1b());
  }
  Cplx result(void) const;
};

} // namespace

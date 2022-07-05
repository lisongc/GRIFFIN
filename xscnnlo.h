/* xscnnlo.h: header file for xscnnlo.cc */

#include "classes.h"

// matel @ NLO is same as matel @ LO
#define mat_SMNLO matel

class mat_SMNNLO : public matel {
public:
  using matel::matel;
  Cplx coeffR(void) const;
  Cplx coeffS1f(void) const;  // 1-loop corrections to S with closed fermion loops
  Cplx coeffS1b(void) const;  // 1-loop corrections without closed fermion loops 
  Cplx coeffS(void) const
  {
    return(matel::coeffS()+coeffS1f()+coeffS1b());
  }
  Cplx result(void) const;
};


class msq_SMNNLO : public matelsq {
public:
  using matelsq::matelsq;
  Cplx result(void) const;
};

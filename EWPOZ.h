/* EWPOZ.h: header file for EWPOZ.cc */

#ifndef __EWPOZ__
#define __EWPOZ__

#include "classes.h"

// axial-vector form factor predicted in the SM (at NLO)
class FA_SMNLO : public FA_SMLO {
public:
  using FA_SMLO::FA_SMLO;
  double res1f(void) const;  // 1-loop corrections with closed fermion loops
  double res1b(void) const;  // 1-loop corrections without closed fermion loops 
  Cplx result(void) const
  {
    return(FA_SMLO::result()+res1f()+res1b());
  }
};

// effective weak mixing angle predicted in the SM (at NLO)
class SW_SMNLO : public SW_SMLO {
public:
  using SW_SMLO::SW_SMLO;
  double res1f(void) const;  // 1-loop corrections with closed fermion loops
  double res1b(void) const;  // 1-loop corrections without closed fermion loops 
  Cplx result(void) const
  {
    return(SW_SMLO::result()+res1f()+res1b());
  }
};

#endif // __EWPOZ__

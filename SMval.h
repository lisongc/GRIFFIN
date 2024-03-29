/* SMval.h: define input class that converts PDG masses to complex-pole masses
   for MW,MZ */

#ifndef __SMval__
#define __SMval__

#include "classes.h"

/* relabel the input parameters:
   - the parameters with indices 0,1,17,18 (which are used in the internal
     computation routines) are now the gauge-boson mass and widths defined in 
     the complex-mass scheme
   - new parameter indices 20-23 are defined for the running-width scheme, 
     which correspond to the masses and width reported in the PDG etc.;
     the symbols MW,MZ,GamW,GamZ are now defined to refer to these */

#define MWc 0
#define MZc 1
#define GWc 17
#define GZc 18

#undef MW
#undef MZ
#undef GamW
#undef GamZ

#define MW 20
#define MZ 21
#define GamW 22
#define GamZ 23

/* input parameter class for SM computations in the complex pole scheme;
   it computes the complex poles masses from user-provided masses in the
   running-width scheme */

class SMval : public inval {
protected:
  void compute(void)
  {
    data[MWc] = data[MW]/sqrt(1+sqr(data[GamW]/data[MW]));
    data[GWc] = data[GamW]/sqrt(1+sqr(data[GamW]/data[MW]));
    data[MZc] = data[MZ]/sqrt(1+sqr(data[GamZ]/data[MZ]));
    data[GZc] = data[GamZ]/sqrt(1+sqr(data[GamZ]/data[MZ]));
  }
public:
  using inval::inval;
  SMval(void) : inval() {};
  SMval(const inval& copyfrom) : inval(copyfrom) {};
};

#endif // __SMval__

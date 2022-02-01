#include "classes.h"

#define MWi 20
#define MZi 21
#define GWi 22
#define GZi 23

/* input parameter class for SM computations in the complex pole scheme;
   it computes the complex poles masses from user-provided masses in the
   running-width scheme */

class SMval : public inval {
  void compute(void)
  {
    data[MW] = data[MWi]/sqrt(1+sqr(data[GWi]/data[MWi]));
    data[GamW] = data[GWi]/sqrt(1+sqr(data[GWi]/data[MWi]));
    data[MZ] = data[MZi]/sqrt(1+sqr(data[GZi]/data[MZi]));
    data[GamZ] = data[GZi]/sqrt(1+sqr(data[GZi]/data[MZi]));
  }
public:
  using inval::inval;
};

#define MWc 0
#define MZc 1
#define GWc 17
#define GZc 18

#undef MW
#undef MZ
#undef GamW
#undef GamZ

#define MW MWi
#define MZ MZi
#define GamW GWi
#define GamZ GZi

/* li.h: header file for classes in li.cc */

#ifndef __li__
#define __li__

#include "Cplx.h"

namespace griffin {

// Li_2: li2(Cplx x), li2(Double x)=li2(x-I*eps), li2c(Double x)=li2(x+I*eps)
Cplx li2(Double x);
Cplx li2c(Double x);
Cplx li2(Cplx x);

} // namespace

#endif

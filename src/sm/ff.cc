/*-----------------------------------------------------------------------------
ff.cc 
Lisong Chen (lic114@pitt.edu), Ayres Freitas (afreitas@pitt.edu)
last revision: 15 Aug 2024
-------------------------------------------------------------------------------
one-loop self-energy, vertex and box form factors
-----------------------------------------------------------------------------*/

#include <math.h>
#include "griffin/sm/ff.h"
#include "griffin/numerics/oneloop.h"
#include "griffin/numerics/li.h"

namespace griffin {

// import automatically generated code
#include "expressions/ff1.in"
#include "expressions/ffs1.in"

extern const double Qf[20];
extern const double I3f[20];

// import automatically generated code
#include "expressions/box1.in"
#include "expressions/boxs1.in"
#include "expressions/box1s0.in"

} // namespace

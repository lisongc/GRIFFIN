/*-----------------------------------------------------------------------------
ff.cc 
Lisong Chen (lic114@pitt.edu), Ayres Freitas (afreitas@pitt.edu)
last revision: 11 Jun 2026
-------------------------------------------------------------------------------
one-loop self-energy, vertex and box form factors
-----------------------------------------------------------------------------*/

#include <math.h>
#include "ff.h"
#include "oneloop.h"
#include "li.h"

namespace griffin {

// import automatically generated code
#include "ff1.in"
#include "ffs1.in"

extern const double Qf[20];
extern const double I3f[20];

// import automatically generated code
#include "box1.in"
#include "boxs1.in"
#include "box1s0.in"

// import SMEFT corrections
#include "ffSMEFT0.in"

} // namespace

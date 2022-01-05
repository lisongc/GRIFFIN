/*-----------------------------------------------------------------------------
ff.cc 
Lisong Chen (lic114@pitt.edu), Ayres Freitas (afreitas@pitt.edu)
last revision: 23 Nov 2021
-------------------------------------------------------------------------------
one-loop self-energy, vertex and box form factors
-----------------------------------------------------------------------------*/

#include <math.h>
#include "ff.h"
#include "oneloop.h"

#include "ff1.in"


//namespace loopt {
//  #include "clooptools.h"
//}
//#define D0 loopt::D0

#include "li.h"

extern const double Qf[5];
extern const double I3f[5];

#include "box1.in"

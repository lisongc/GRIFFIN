/*-----------------------------------------------------------------------------
ff.cc 
Lisong Chen (lic114@pitt.edu), Ayres Freitas (afreitas@pitt.edu)
last revision: 18 Oct 2022
-------------------------------------------------------------------------------
one-loop self-energy, vertex and box form factors
-----------------------------------------------------------------------------*/

#include <math.h>
#include "ff.h"
#include "oneloop.h"

// import automatically generated code
#include "ff1.in"
#include "ffs1.in"

#include "li.h"

extern const double Qf[20];
extern const double I3f[20];

// import automatically generated code
#include "box1.in"
#include "boxs1.in"

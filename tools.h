/* tools.h: header file for tools.cc */

#include "classes.h"

#define COMPPOLESCHEME  0
#define RUNWIDTHSCHEME  1

// Compute partial width for Z->ff using form factor objects provided, and
//  additionally including final-state QCD/QED corrections
// Note: partzwidth will modify the FA and FV objects passed to it,
//        by setting the fermion type and input
// type: fermion type f
// scheme: definition of the output value for the width,
//          can take values COMPPOLESCHEME (for complex-pole scheme)
//          or RUNWIDTHSCHEME (for the running-width scheme)
double partzwidth(FA_SMLO& fa, FV_SMLO& fv, const int type, 
                  const inval& input, int scheme = COMPPOLESCHEME);

// Compute total width for Z->ff using form factor objects provided
// Note: zwidth will modify the FA and FV objects passed to it,
//        by setting the fermion type and input
double zwidth(FA_SMLO& fa, FV_SMLO& fv, const inval& input,
	      int scheme = COMPPOLESCHEME);

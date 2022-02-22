/* tools.h: header file for tools.cc */

#include "classes.h"

#define COMPPOLESCHEME  0
#define RUNWIDTHSCHEME  1

// Compute partial width for Z->ff using form factors provided
// Note: partzwidth will modify the FA and FV objects passed to it,
//        by setting the fermion type and input
double partzwidth(FA_SMLO& fa, FV_SMLO& fv, const int type, const int gen,
                  const inval& input, int scheme = COMPPOLESCHEME);

// Compute partial width for Z->ff using form factors provided
// Note: zwidth will modify the FA and FV objects passed to it,
//        by setting the fermion type and input
double zwidth(FA_SMLO& fa, FV_SMLO& fv, const inval& input,
	      int scheme = COMPPOLESCHEME);

/*
  File: delrho.h
  Purpose: Declare higher-order rho-parameter and running strong-coupling functions.
  Typed-input integration: Jackson Wallace and Lisong Chen
  Last updated: 2026-11-16
  Original GRIFFIN calculations/author records are retained below or in the
  including source file; this migration changes labels, not the expressions.
*/

#ifndef GRIFFIN_SM_DELRHO_H
#define GRIFFIN_SM_DELRHO_H

/* delrho.h: header file for delrho.cc */

// alpha scheme for leading mt-enhanced contributions:
#define DRHOSCHEME  1

// Gmu scheme for leading mt-enhanced contributions:
//#define DRHOSCHEME  (0.4501581580785531*ival->get(InputPar::fermiConstant)*MWs*(1-MWs/MZs)/ival->get(InputPar::alpha))

#include "griffin/core/classes.h"

namespace griffin {

double getalphas(double musqval, double mzsqval, double alphasmz);

double delrho2a2(const inval* ival);   // yt^4 corrections to \Delta\rho
double delrho2aas(const inval* ival);  // yt^2*as corrections to \Delta\rho
double delrho3a3(const inval* ival);   // yt^6 corrections to \Delta\rho
double delrho3a2as(const inval* ival); // yt^4*as corrections to \Delta\rho
double delrho3aas2(const inval* ival); // yt^2*as^2 corrections to \Delta\rho
double delrho4aas3(const inval* ival); // yt^2*as^3 corrections to \Delta\rho

}

#endif // GRIFFIN_SM_DELRHO_H

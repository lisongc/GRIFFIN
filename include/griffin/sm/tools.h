/*
  File: tools.h
  Purpose: Declare decay-width utilities and conversion of supplied amplitudes
           into a differential cross section.
  Original utility declarations: GRIFFIN contributors; see src/sm/tools.cc.
  Integration contributors: Jackson Wallace and Lisong Chen
  Last updated: 2026-11-16
*/

#ifndef GRIFFIN_SM_TOOLS_H
#define GRIFFIN_SM_TOOLS_H

/* tools.h: header file for tools.cc */

#include <array>

#include "griffin/core/classes.h"

namespace griffin {

#define COMPPOLESCHEME  0
#define RUNWIDTHSCHEME  1

// Convert four supplied reduced amplitudes, ordered VV, AV, VA, AA, into
// unpolarized d sigma / d cos(theta) for massless external fermions.
// Supports e+e- -> mu+mu-, tau+tau- and u,d,s,c,b quark pairs; includes the
// final-state color factor (1 or 3). Bhabha and neutrino channels are rejected.
// This is only the angular square: it neither constructs amplitudes nor adds
// SM/SMEFT corrections. The input amplitudes must use GRIFFIN's normalization.
// s is in GeV^2; output is GeV^-2 by default (unitConversion=0.38937966e9 for pb).
// Requires finite s>0, -1<=cosTheta<=1 and unitConversion>0. Whether the
// amplitudes themselves are valid at a given energy/angle is the caller's
// responsibility; this utility does not repair near-pole/collinear limits.
double differentialCrossSectionFromAmplitudes(
    const std::array<Cplx,4>& amplitudes, Fermion initial, Fermion final,
    double invariantMassSquared, double cosTheta, double unitConversion = 1.0);

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

} // namespace

#endif // GRIFFIN_SM_TOOLS_H

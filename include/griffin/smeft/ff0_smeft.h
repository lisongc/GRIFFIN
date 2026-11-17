/*
  File: ff0_smeft.h
  Purpose: Declare direct SMEFT changes to the axial/vector Zff couplings,
           selected only by typed Fermion labels.
  Authors: Jackson Wallace and Lisong Chen
  Last updated: 2026-11-16

  These functions return the explicit dimension-six tree vertex shifts, not
  the SM couplings themselves.  The supplied alpha, MW and MZ are held fixed
  while evaluating the formula.  There is no input-scheme selector and no
  additional W-mass response in these helpers.

  SMEFTCalculator supplies its already-resolved input point.  In Gmu mode
  that point contains the finite coefficient-dependent W mass.  Keeping the
  direct formula linear at a fixed point does not linearize the subsequent
  pole construction or the cross section.
*/

#ifndef GRIFFIN_SMEFT_FF0_H
#define GRIFFIN_SMEFT_FF0_H

#include "griffin/core/classes.h"
#include "griffin/smeft/SMEFTInput.h"

namespace griffin {
namespace smeft {

// Explicit change in the axial Zff coupling at the supplied mass point.
double directZAxialShift(Fermion type, const inval& smInput,
                        const Input& eftInput);

// Explicit change in the vector Zff coupling at the supplied mass point.
double directZVectorShift(Fermion type, const inval& smInput,
                         const Input& eftInput);

} // namespace smeft
} // namespace griffin

#endif // GRIFFIN_SMEFT_FF0_H

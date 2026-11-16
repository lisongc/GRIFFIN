/*
  File: SMEFT.h
  Purpose: Collect the finite SMEFT interface and SM input headers in a safe order.
  Authors: Jackson Wallace and Lisong Chen
  Last updated: 2026-11-16
*/

#ifndef GRIFFIN_SMEFT_H
#define GRIFFIN_SMEFT_H

// Recommended finite-SMEFT interface and SM input containers.
// Parse SM physics headers before the input containers: legacy SM input
// headers redefine short mass macros. Include this entry header first when
// mixing it with advanced SM headers; use explicit InputPar names in callers.
#include "griffin/sm/EWPOZ2.h"
#include "griffin/smeft/SMEFTCalculator.h"
#include "griffin/sm/SMvalG.h"

#endif // GRIFFIN_SMEFT_H

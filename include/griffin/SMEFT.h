#ifndef GRIFFIN_SMEFT_H
#define GRIFFIN_SMEFT_H

// Recommended finite-SMEFT interface and SM input containers.
// Parse SM physics headers before the input containers: legacy SM input
// headers redefine short mass macros. Include this entry header first when
// mixing it with advanced SM headers; use explicit InputPar names in callers.
#include "griffin/sm/EWPOZ2.h"
#include "griffin/SMEFTCalculator.h"
#include "griffin/sm/SMvalG.h"

#endif // GRIFFIN_SMEFT_H

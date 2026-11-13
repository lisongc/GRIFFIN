/*
  NativeInputPoint.h: finite-coefficient electroweak input resolution.

  The Wilson container remains independent of the electroweak inputs.  This
  object owns a plain inval copy, so updating its derived W parameters cannot
  invoke SMvalGmu::compute() or mutate the caller's input object. This is the
  mass-resolution record used by SMEFTCalculator, not an alternative evaluator.
*/

#ifndef GRIFFIN_SMEFT_NATIVE_INPUT_POINT_H
#define GRIFFIN_SMEFT_NATIVE_INPUT_POINT_H

#include "griffin/core/classes.h"
#include "griffin/EWInputScheme.h"
#include "griffin/smeft/SMEFTInput.h"

namespace griffin {
namespace smeft {

struct NativeSolverOptions {
  // Stopping criterion for the residual of the equation for MW^2, not MW.
  double massSquaredToleranceGeV2 = 1.0e-4;
  unsigned maxIterations = 1000;
};

class NativeInputPoint {
public:
  /*
    The supplied pole Z parameters are authoritative and are not reconverted.
    AlphaMWMZ returns an unchanged copy of the supplied SM input.  GmuAlphaMZ
    solves Jackson's finite fixed-point equation, with the dimension-six
    additive shift outside the SM square root, and updates all four W slots.

    At convergence the stored mass is the last mass at which Delta r_SM was
    evaluated, not the subsequent candidate.  This deliberately matches the
    original Jackson solver's endpoint convention at the default tolerance.
  */
  NativeInputPoint(const inval& smInput, const Input& eftInput,
                   EWInputScheme scheme,
                   const NativeSolverOptions& options = NativeSolverOptions());

  const inval& smInput() const { return smInput_; }
  // Dimensionless additive term in the MW^2/MZ^2 fixed-point equation.
  // This is not the final difference MW(C)^2 - MW(0)^2 divided by MZ^2.
  double additiveMassSquaredShiftOverMZ2() const { return additiveShift_; }
  unsigned iterations() const { return iterations_; }

  // |fixedPointMap(MW^2) - MW^2| at the returned point; zero for AlphaMWMZ.
  double finalResidualGeV2() const { return finalResidualGeV2_; }

private:
  inval smInput_;
  double additiveShift_;
  unsigned iterations_;
  double finalResidualGeV2_;
};

} // namespace smeft
} // namespace griffin

#endif // GRIFFIN_SMEFT_NATIVE_INPUT_POINT_H

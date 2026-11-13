/*
  WMassSolution.h: the W-mass solution and the electroweak inputs that use it.

  In GmuAlphaMZ, MW is a prediction rather than an independent input. This
  object solves Jackson's finite mass equation using the supplied SM inputs
  and Wilson coefficients, then updates the associated W width and converted
  W mass/width values. It does not use a linearized W-mass response.
  In AlphaMWMZ, MW is already supplied: the input copy is checked and retained
  without solving a mass equation.

  The object owns a plain inval copy. Updating its W parameters therefore
  neither changes the caller's input nor triggers SMvalGmu::compute(), which
  could otherwise replace the SMEFT mass solution with an SM-only solution.
  It also records the iteration count and final mass-equation residual.

  SMEFTCalculator constructs two solutions automatically: one with the user's
  coefficients, and one with zero coefficients for the SM reference. Ordinary
  users need not construct this class or include this header directly; use
  <griffin/SMEFT.h>. Optional solver diagnostics are available through
  calculator.massSolution(). This is not another prediction engine or an
  energy/angle point for scattering. The solving code is in WMassSolution.cc.
*/

#ifndef GRIFFIN_SMEFT_W_MASS_SOLUTION_H
#define GRIFFIN_SMEFT_W_MASS_SOLUTION_H

#include "griffin/core/classes.h"
#include "griffin/EWInputScheme.h"
#include "griffin/smeft/SMEFTInput.h"

namespace griffin {
namespace smeft {

struct WMassSolverOptions {
  // Stopping criterion for the residual of the equation for MW^2, not MW.
  double massSquaredToleranceGeV2 = 1.0e-4;
  unsigned maxIterations = 1000;
};

class WMassSolution {
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
  WMassSolution(const inval& smInput, const Input& eftInput,
                EWInputScheme scheme,
                const WMassSolverOptions& options = WMassSolverOptions());

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

#endif // GRIFFIN_SMEFT_W_MASS_SOLUTION_H

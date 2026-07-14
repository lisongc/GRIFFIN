/*
  EWInputScheme.h: electroweak input-scheme support for SMEFT calculations

  The Wilson coefficients and Lambda stored by smeft::Input are independent of
  the electroweak input scheme.  This module describes how the accompanying SM
  input is interpreted when those coefficients are converted into physical
  dimension-six corrections.

  Two schemes are currently supported:

    GmuAlphaMZ  inputs {G_mu, alpha, M_Z}; M_W is the SM value derived by
                  SMvalGmu.  This is often called the LEP scheme.

    AlphaMWMZ   inputs {alpha, M_W, M_Z}; M_W is supplied directly through
                  SMval, while G_mu (or Delta r) is a prediction.

  Explicit parameter-set names are used because the shorter scheme names are
  not uniform between the SM and SMEFT literature.

  The AlphaMWMZ definitions follow the tree-level relations in Eqs. (2.8) and
  (2.11) of arXiv:2305.03763, adapted to this package's convention of
  dimensionless C_i with an explicit 1/Lambda^2.  In GmuAlphaMZ the stored MW
  is the SM solution and its dimension-six displacement is represented by the
  explicitly linear DerivedMWShiftD6 below; the SM input object is never
  mutated by Wilson coefficients.
*/

#ifndef GRIFFIN_SMEFT_EW_INPUT_SCHEME_H
#define GRIFFIN_SMEFT_EW_INPUT_SCHEME_H

#include "classes.h"
#include "prediction/EWInputScheme.h"
#include "smeft/SMEFTInput.h"

namespace griffin {
namespace smeft {

// Compatibility alias.  The input scheme now belongs to the shared
// prediction context, while existing smeft::EWInputScheme source continues to
// compile unchanged.
using ::griffin::EWInputScheme;

// Stable human-readable name for diagnostics and test output.
const char* inputSchemeName(EWInputScheme scheme);

/*
  Linear dimension-six displacement of the MW value derived from
  {G_mu, alpha, M_Z}.  inputDeltaMW2OverMZ2 is the additive quantity called
  resSMEFT() in Jackson's fixed-point equation.  Since the SM Delta r in that
  equation also changes when MW changes, responseFactor converts this input
  shift into the final derived-mass displacement deltaMW2OverMZ2.

  linearizedMW() is a convenient prediction, but it must not be fed back into
  an observable that already uses the GmuAlphaMZ vertex shifts: those shifts
  contain the same first-order mass displacement.
*/
struct DerivedMWShiftD6 {
  double referenceMW;
  double inputDeltaMW2OverMZ2;
  double responseFactor;
  double deltaMW2OverMZ2;
  double deltaMW2;
  double deltaMW;

  double linearizedMW() const { return referenceMW+deltaMW; }
};

DerivedMWShiftD6 derivedMWShiftD6(const inval& smInput,
                                  const Input& eftInput);

namespace detail {

/*
  Leading electroweak parameters used to assemble a linear dimension-six term.

  The suffix 0 emphasizes that e0 and gZ0 are reference SM couplings.  SMEFT
  corrections must not be absorbed into these reference values.  epsilonD6 is
  the scheme reference v^2 divided by Lambda^2; replacing the exact SMEFT vev
  by that reference is valid through O(1/Lambda^2).
*/
struct EWTreeParameters {
  double e0;
  double sw0;
  double cw0;
  double gZ0;
  double vevSquared;
  double epsilonD6;
};

/*
  Construct the common tree-level context.

    GmuAlphaMZ: v^2 = 1/(sqrt(2) G_mu)
    AlphaMWMZ:  v^2 = (2 M_W s_W/e)^2

  In both current schemes c_W=M_W/M_Z and e=sqrt(4 pi alpha).  The origin of
  M_W differs: it is derived in GmuAlphaMZ and directly supplied in AlphaMWMZ.
*/
EWTreeParameters ewTreeParameters(EWInputScheme scheme,
                                  const inval& smInput,
                                  const Input& eftInput);

/*
  Muon-decay Wilson combination

    Delta_mu = Cphil3[1,1] + Cphil3[2,2]
             - (Cll[1,2,2,1] + Cll[2,1,1,2])/2.

  It is shared by the GmuAlphaMZ input transformation and the AlphaMWMZ
  prediction for Delta r.  Keeping one implementation prevents the flavor
  convention from drifting between the two calculations.
*/
double deltaMu(const Input& input);

// Dimension-six tree contribution to Delta r in the selected input scheme.
// This common implementation is also used to derive delta(MW^2/MZ^2), so the
// two representations cannot drift apart coefficient by coefficient.
double deltaRSMEFTTree(EWInputScheme scheme,
                       const inval& smInput,
                       const Input& eftInput);

} // namespace detail
} // namespace smeft
} // namespace griffin

#endif // GRIFFIN_SMEFT_EW_INPUT_SCHEME_H

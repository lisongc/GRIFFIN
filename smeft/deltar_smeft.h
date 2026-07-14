/*
  deltar_smeft.h: Delta r at SM NLO plus SMEFT tree level

  Delta r relates the two reference vevs

    v_alpha^2 / v_mu^2 = 1 + Delta r,

  where v_mu^2=1/(sqrt(2) G_mu) and
  v_alpha=2 M_W s_W/sqrt(4 pi alpha).

  The interpretation depends on EWInputScheme:

    AlphaMWMZ: M_W is an input and Delta r (equivalently G_mu) is a
                prediction.  This is the natural scheme for using result()
                as an EW pseudo-observable.

    GmuAlphaMZ: G_mu is an input.  The raw SMEFT input-transformation
                coefficient uses the tree-level relation
                  sw^2 cw^2 = pi alpha/(sqrt(2) G_mu M_Z^2)
                and derivedMWShiftD6() propagates it through the same
                loop-improved SM fixed point used for the stored M_W.  It may
                be inspected, but it is not independent of the derived mass.

  Following the cumulative convention of the other SMEFT pseudo-observables,

    dr_SMEFTLO::result() = dr_SMNLO::result() + resD6Tree().

  No SMEFT-shifted M_W is written into the associated inval.  In particular,
  callers must not feed an independently shifted M_W back into LEP-scheme Z
  form factors that already contain the corresponding universal correction.
*/

#ifndef GRIFFIN_SMEFT_DELTAR_H
#define GRIFFIN_SMEFT_DELTAR_H

#include "deltar.h"
#include "smeft/EWInputScheme.h"
#include "smeft/SMEFTInput.h"

namespace griffin {
namespace smeft {

class dr_SMEFTLO : public dr_SMNLO {
private:
  const Input* eftInput;
  EWInputScheme scheme;

public:
  /*
    The scheme argument is mandatory because the same numerical Delta-r
    relation has a different input/output role in the two schemes.
  */
  dr_SMEFTLO(const inval& smInput, const Input& input,
             EWInputScheme inputScheme)
    : dr_SMNLO(smInput), eftInput(&input), scheme(inputScheme) {}

  void setSMEFTInput(const Input& input) { eftInput = &input; }
  EWInputScheme inputScheme() const { return scheme; }

  // Pure O(1/Lambda^2) tree-level contribution.
  Cplx resD6Tree() const;
  // Full SM NLO plus the linear dimension-six tree contribution.
  Cplx result() const override;
};

} // namespace smeft
} // namespace griffin

#endif // GRIFFIN_SMEFT_DELTAR_H

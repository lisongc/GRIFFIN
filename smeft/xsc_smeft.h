/*
  xsc_smeft.h: matrix elements at SM LO plus linear SMEFT tree level

  mat_SMEFTLO follows GRIFFIN's cumulative-order naming convention:

    result() = SM LO + terms linear in dimension-six Wilson coefficients.

  The supplied FA and SW objects must contain an SM prediction only.  The
  class adds the Z-vertex shifts explicitly and linearizes the two-vertex
  residue, thereby omitting the unwanted product delta(z_i)*delta(z_f).

  The numerical amplitude is organized as an improved pole plus an exact
  regular background:

    resultD6Tree() = poleD6Tree() + resoffZD6Tree().

  The S and S' methods expose Laurent coefficients for pole fits and checks;
  result() never adds them separately.  Local chirality-preserving
  four-fermion operators obey resoffZ4fD6Tree() = coeffS4fD6Tree() and
  coeffSp4fD6Tree() = 0 because their reduced amplitude is constant.

  The electron-neutrino contact coefficients are available, but the complete
  background and result are deliberately disabled until the physical
  t-channel W contribution is implemented in a dedicated process module.
*/

#ifndef GRIFFIN_SMEFT_XSC_H
#define GRIFFIN_SMEFT_XSC_H

#include "classes.h"
#include "smeft/EWInputScheme.h"
#include "smeft/SMEFTInput.h"

namespace griffin {
namespace smeft {

/*
  Standalone local-contact coefficients.  These functions implement the
  massless vector/axial formulas documented in SMEFT_IMPLEMENTATION.md.
  They assume an electron initial state and reject Bhabha, scalar, and
  pseudoscalar channels, which require separate process structures.
*/
Cplx coeffS4fD6Tree(int intype, int outtype, int inform, int outform,
                    const Input& eftInput);
Cplx coeffSp4fD6Tree(int intype, int outtype, int inform, int outform,
                     const Input& eftInput);

inline Cplx coeffS4fD6Tree(Fermion intype, Fermion outtype, Current inform,
                           Current outform, const Input& eftInput)
{
  return coeffS4fD6Tree(index(intype),index(outtype),index(inform),
                        index(outform),eftInput);
}

inline Cplx coeffSp4fD6Tree(Fermion intype, Fermion outtype, Current inform,
                            Current outform, const Input& eftInput)
{
  return coeffSp4fD6Tree(index(intype),index(outtype),index(inform),
                         index(outform),eftInput);
}

/*
  Observable-level massless differential cross section, truncated exactly at
  O(1/Lambda^2).  sm is the prediction made from the supplied SM-only form
  factors.  d6 is the interference of the tree SM amplitude with the linear
  SMEFT tree amplitude; |M_D6|^2 is never formed.  Values are in GeV^-2 unless
  the caller supplies a unitConversion (for example GeVtoNB).
*/
struct DifferentialCrossSectionLinear {
  double sm;
  double d6;

  double total() const { return sm+d6; }
};

DifferentialCrossSectionLinear differentialCrossSectionLinear(
    int intype, int outtype,
    const psobs& FAin, const psobs& FAout,
    const psobs& SWin, const psobs& SWout,
    double sval, double costheta, const inval& smInput,
    const Input& eftInput,
    EWInputScheme inputScheme = EWInputScheme::GmuAlphaMZ,
    double unitConversion = 1.0);

inline DifferentialCrossSectionLinear differentialCrossSectionLinear(
    Fermion intype, Fermion outtype,
    const psobs& FAin, const psobs& FAout,
    const psobs& SWin, const psobs& SWout,
    double sval, double costheta, const inval& smInput,
    const Input& eftInput,
    EWInputScheme inputScheme = EWInputScheme::GmuAlphaMZ,
    double unitConversion = 1.0)
{
  return differentialCrossSectionLinear(
      index(intype),index(outtype),FAin,FAout,SWin,SWout,sval,costheta,
      smInput,eftInput,inputScheme,unitConversion);
}

class mat_SMEFTLO : public matel {
private:
  const Input* eftInput;
  EWInputScheme scheme;

public:
  // Legacy macro API.  FA/SW objects must be SM-only predictions.
  mat_SMEFTLO(int intype, int outtype, int inform, int outform,
              const psobs& FAin, const psobs& FAout,
              const psobs& SWin, const psobs& SWout,
              double sval, double costheta, const inval& smInput,
              const Input& input,
              EWInputScheme inputScheme = EWInputScheme::GmuAlphaMZ)
    : matel(intype,outtype,inform,outform,FAin,FAout,SWin,SWout,
            sval,costheta,smInput),
      eftInput(&input), scheme(inputScheme) {}

  // Strongly typed API with the same behavior as the legacy constructor.
  mat_SMEFTLO(Fermion intype, Fermion outtype,
              Current inform, Current outform,
              const psobs& FAin, const psobs& FAout,
              const psobs& SWin, const psobs& SWout,
              double sval, double costheta, const inval& smInput,
              const Input& input,
              EWInputScheme inputScheme = EWInputScheme::GmuAlphaMZ)
    : mat_SMEFTLO(index(intype),index(outtype),index(inform),index(outform),
                  FAin,FAout,SWin,SWout,sval,costheta,smInput,input,
                  inputScheme) {}

  void setSMEFTInput(const Input& input) { eftInput = &input; }
  EWInputScheme inputScheme() const { return scheme; }

  // Z-pole residue: SM value plus its explicitly linearized D6 correction.
  Cplx coeffRD6Tree() const;
  Cplx coeffR() const;

  // Local four-fermion contribution to S and its vanishing tree derivative.
  // These are coefficient diagnostics and are not added by result().
  Cplx coeffS4fD6Tree() const;
  Cplx coeffSp4fD6Tree() const;
  Cplx coeffS() const;
  Cplx coeffSp() const;

  /*
    Numerical amplitude decomposition.  resoffZ4fD6Tree() is the exact local
    contact amplitude.  resoffZD6Tree() is the complete regular D6 remainder
    available for this class and is the extension point for future nonlocal
    process-specific contributions.  No S coefficient is added here.
  */
  Cplx poleD6Tree() const;
  Cplx resoffZ4fD6Tree() const;
  Cplx resoffZD6Tree() const;
  Cplx resultD6Tree() const;

  // Cumulative amplitude interface.  Squaring result() directly also forms an
  // O(Lambda^-4) term; use differentialCrossSectionLinear() for a cross section
  // that is structurally truncated at linear dimension-six order.
  Cplx resoffZ() const;
  Cplx result() const override;
};

} // namespace smeft
} // namespace griffin

#endif // GRIFFIN_SMEFT_XSC_H

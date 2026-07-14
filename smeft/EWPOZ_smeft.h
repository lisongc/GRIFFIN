/*
  EWPOZ_smeft.h: Z-pole pseudo-observables at SM LO + SMEFT tree level

  The class names follow GRIFFIN's cumulative-order convention:

    *_SMEFTLO::result() = SM LO + linear dimension-six tree contribution.

  resD6Tree() exposes the correction by itself for diagnostics, fits, and
  future combinations with higher-order SM predictions.  These classes do not
  modify the pole-residue or regular matrix-element coefficients.

  Each object stores an EWInputScheme.  The SM input paired with it must obey
  the corresponding contract:

    GmuAlphaMZ -> use SMvalGmu, whose stored MW is the SM-derived value;
    AlphaMWMZ  -> use SMval, with MW supplied directly.

  The old three-argument constructors remain GmuAlphaMZ for compatibility.
*/

#ifndef GRIFFIN_SMEFT_EWPOZ_H
#define GRIFFIN_SMEFT_EWPOZ_H

#include "classes.h"
#include "smeft/EWInputScheme.h"
#include "smeft/SMEFTInput.h"

namespace griffin {
namespace smeft {

class SW_SMEFTLO : public SW_SMLO {
private:
  const Input* eftInput;
  EWInputScheme scheme;

public:
  SW_SMEFTLO(int type, const inval& smInput, const Input& input,
             EWInputScheme inputScheme = EWInputScheme::GmuAlphaMZ)
    : SW_SMLO(type,smInput), eftInput(&input), scheme(inputScheme) {}

  SW_SMEFTLO(Fermion type, const inval& smInput, const Input& input,
             EWInputScheme inputScheme = EWInputScheme::GmuAlphaMZ)
    : SW_SMEFTLO(index(type),smInput,input,inputScheme) {}

  void setSMEFTInput(const Input& input) { eftInput = &input; }
  EWInputScheme inputScheme() const { return scheme; }
  Cplx resD6Tree() const;
  Cplx result() const override;
};

class FA_SMEFTLO : public FA_SMLO {
private:
  const Input* eftInput;
  EWInputScheme scheme;

public:
  FA_SMEFTLO(int type, const inval& smInput, const Input& input,
             EWInputScheme inputScheme = EWInputScheme::GmuAlphaMZ)
    : FA_SMLO(type,smInput), eftInput(&input), scheme(inputScheme) {}
  FA_SMEFTLO(Fermion type, const inval& smInput, const Input& input,
             EWInputScheme inputScheme = EWInputScheme::GmuAlphaMZ)
    : FA_SMEFTLO(index(type),smInput,input,inputScheme) {}

  void setSMEFTInput(const Input& input) { eftInput = &input; }
  EWInputScheme inputScheme() const { return scheme; }
  Cplx resD6Tree() const;
  Cplx result() const override;
};

class FV_SMEFTLO : public FV_SMLO {
private:
  const Input* eftInput;
  EWInputScheme scheme;

public:
  FV_SMEFTLO(int type, const inval& smInput, const Input& input,
             EWInputScheme inputScheme = EWInputScheme::GmuAlphaMZ)
    : FV_SMLO(type,smInput), eftInput(&input), scheme(inputScheme) {}

  FV_SMEFTLO(Fermion type, const inval& smInput, const Input& input,
             EWInputScheme inputScheme = EWInputScheme::GmuAlphaMZ)
    : FV_SMEFTLO(index(type),smInput,input,inputScheme) {}

  void setSMEFTInput(const Input& input) { eftInput = &input; }
  EWInputScheme inputScheme() const { return scheme; }
  Cplx resD6Tree() const;
  Cplx result() const override;
};

} // namespace smeft
} // namespace griffin

#endif // GRIFFIN_SMEFT_EWPOZ_H

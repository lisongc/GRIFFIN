/*
  ElectroweakSetup.h
  Purpose: define the electroweak input-scheme choices and pair a chosen
  scheme with compatible ordinary SM inputs. This descriptor does not
  calculate observables or contain Wilson coefficients.
  Merged-version authors: Jackson Wallace and Lisong Chen.
  Last updated: 2026-11-16.
*/

#ifndef GRIFFIN_SM_ELECTROWEAK_SETUP_H
#define GRIFFIN_SM_ELECTROWEAK_SETUP_H

#include <stdexcept>

#include "griffin/core/classes.h"
#include "griffin/sm/deltar.h"

namespace griffin {

// Shared electroweak input-scheme identifiers. Use griffin::EWInputScheme;
// there is no separate SMEFT-namespace alias.
enum class EWInputScheme {
  GmuAlphaMZ,
  AlphaMWMZ
};

inline const char* ewInputSchemeName(EWInputScheme scheme)
{
  switch(scheme)
  {
    case EWInputScheme::GmuAlphaMZ:
      return "{Gmu, alpha, MZ}";
    case EWInputScheme::AlphaMWMZ:
      return "{alpha, MW, MZ}";
  }
  return "unknown electroweak input scheme";
}

namespace prediction {

// The electroweak input scheme paired with its ordinary SM inputs.
// This lightweight descriptor borrows its input, which must remain alive
// while the descriptor is used. SMEFTCalculator copies that input at
// construction and does not retain the caller's setup or input pointer.
class ElectroweakSetup {
private:
  const inval* input_;
  EWInputScheme scheme_;

public:
  ElectroweakSetup(EWInputScheme scheme, const inval& input)
    : input_(&input), scheme_(scheme)
  {
    const bool derivesMWFromGmu =
        dynamic_cast<const invalGmu*>(&input) != 0;
    if(scheme == EWInputScheme::GmuAlphaMZ && !derivesMWFromGmu)
      throw std::invalid_argument(
          "GmuAlphaMZ requires invalGmu/SMvalGmu SM input values");
    if(scheme == EWInputScheme::AlphaMWMZ && derivesMWFromGmu)
      throw std::invalid_argument(
          "AlphaMWMZ requires SM input values with a directly supplied MW");
  }

  const inval& input() const { return *input_; }
  EWInputScheme inputScheme() const { return scheme_; }
};

} // namespace prediction
} // namespace griffin

#endif // GRIFFIN_SM_ELECTROWEAK_SETUP_H

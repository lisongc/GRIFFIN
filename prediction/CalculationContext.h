#ifndef GRIFFIN_PREDICTION_CALCULATION_CONTEXT_H
#define GRIFFIN_PREDICTION_CALCULATION_CONTEXT_H

#include <stdexcept>

#include "classes.h"
#include "deltar.h"
#include "prediction/EWInputScheme.h"

namespace griffin {
namespace prediction {

// A single context is shared by the SM reference and every model provider.
// Pointer identity is intentional: it prevents two independently mutable
// input objects from being mistaken for the same reference point.
class CalculationContext {
private:
  const inval* input_;
  EWInputScheme scheme_;

public:
  CalculationContext(EWInputScheme scheme, const inval& input)
    : input_(&input), scheme_(scheme)
  {
    const bool derivesMWFromGmu =
        dynamic_cast<const invalGmu*>(&input) != 0;
    if(scheme == EWInputScheme::GmuAlphaMZ && !derivesMWFromGmu)
      throw std::invalid_argument(
          "GmuAlphaMZ requires an invalGmu/SMvalGmu reference input");
    if(scheme == EWInputScheme::AlphaMWMZ && derivesMWFromGmu)
      throw std::invalid_argument(
          "AlphaMWMZ requires a directly supplied MW reference input");
  }

  const inval& input() const { return *input_; }
  const inval* inputAddress() const { return input_; }
  EWInputScheme inputScheme() const { return scheme_; }

  bool compatibleWith(const CalculationContext& other) const
  {
    return input_ == other.input_ && scheme_ == other.scheme_;
  }
};

} // namespace prediction
} // namespace griffin

#endif // GRIFFIN_PREDICTION_CALCULATION_CONTEXT_H

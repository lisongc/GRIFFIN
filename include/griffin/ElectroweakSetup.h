#ifndef GRIFFIN_PREDICTION_ELECTROWEAK_SETUP_H
#define GRIFFIN_PREDICTION_ELECTROWEAK_SETUP_H

#include <stdexcept>

#include "griffin/core/classes.h"
#include "griffin/sm/deltar.h"
#include "griffin/EWInputScheme.h"

namespace griffin {
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
          "GmuAlphaMZ requires an invalGmu/SMvalGmu reference input");
    if(scheme == EWInputScheme::AlphaMWMZ && derivesMWFromGmu)
      throw std::invalid_argument(
          "AlphaMWMZ requires a directly supplied MW reference input");
  }

  const inval& input() const { return *input_; }
  EWInputScheme inputScheme() const { return scheme_; }
};

} // namespace prediction
} // namespace griffin

#endif // GRIFFIN_PREDICTION_ELECTROWEAK_SETUP_H

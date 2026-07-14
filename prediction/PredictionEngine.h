#ifndef GRIFFIN_PREDICTION_ENGINE_H
#define GRIFFIN_PREDICTION_ENGINE_H

#include "prediction/CombinationRule.h"
#include "prediction/PredictionResult.h"
#include "prediction/TheoryProvider.h"

namespace griffin {
namespace prediction {

class PredictionEngine {
private:
  const TheoryProvider* reference_;
  const TheoryProvider* model_;
  CombinationRule combination_;

  void validateContexts() const;

public:
  explicit PredictionEngine(const TheoryProvider& reference)
    : reference_(&reference), model_(0),
      combination_(CombinationRule::additiveCorrection())
  {
    validateContexts();
  }

  PredictionEngine(const TheoryProvider& reference,
                   const TheoryProvider& model,
                   CombinationRule combination)
    : reference_(&reference), model_(&model), combination_(combination)
  {
    validateContexts();
  }

  const CalculationContext& context() const { return reference_->context(); }
  const CombinationRule& combinationRule() const { return combination_; }

  PredictionResult<Cplx> pseudoObservable(
      PseudoObservable observable, Fermion fermion,
      const Accuracy& accuracy) const;

  PredictionResult<double> differentialCrossSection(
      const ScatteringPoint& point,
      const Accuracy& accuracy) const;
};

} // namespace prediction
} // namespace griffin

#endif // GRIFFIN_PREDICTION_ENGINE_H

#ifndef GRIFFIN_PREDICTION_RESULT_H
#define GRIFFIN_PREDICTION_RESULT_H

#include "prediction/PerturbativeSeries.h"

namespace griffin {
namespace prediction {

template<typename Value>
class PredictionResult {
private:
  Value reference_;
  Value modelCorrection_;
  PerturbativeSeries<Value> referenceSeries_;
  PerturbativeSeries<Value> modelSeries_;

public:
  PredictionResult(const Value& reference, const Value& modelCorrection,
                   const PerturbativeSeries<Value>& referenceSeries,
                   const PerturbativeSeries<Value>& modelSeries)
    : reference_(reference), modelCorrection_(modelCorrection),
      referenceSeries_(referenceSeries), modelSeries_(modelSeries) {}

  const Value& reference() const { return reference_; }
  const Value& modelCorrection() const { return modelCorrection_; }
  Value total() const { return reference_+modelCorrection_; }

  const PerturbativeSeries<Value>& referenceSeries() const
  {
    return referenceSeries_;
  }

  const PerturbativeSeries<Value>& modelSeries() const
  {
    return modelSeries_;
  }
};

} // namespace prediction
} // namespace griffin

#endif // GRIFFIN_PREDICTION_RESULT_H

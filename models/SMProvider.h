#ifndef GRIFFIN_MODELS_SM_PROVIDER_H
#define GRIFFIN_MODELS_SM_PROVIDER_H

#include "prediction/TheoryProvider.h"

namespace griffin {
namespace prediction {

class SMProvider : public TheoryProvider {
private:
  CalculationContext context_;

public:
  explicit SMProvider(const CalculationContext& context)
    : context_(context) {}

  const CalculationContext& context() const override { return context_; }
  ProviderOutputKind outputKind() const override
  {
    return ProviderOutputKind::ReferenceTheory;
  }
  const char* name() const override { return "Standard Model"; }

  PerturbativeSeries<Cplx> pseudoObservableSeries(
      PseudoObservable observable, Fermion fermion) const override;

  PerturbativeSeries<Cplx> amplitudeSeries(
      const ScatteringPoint& point,
      const CurrentChannel& channel) const override;
};

} // namespace prediction
} // namespace griffin

#endif // GRIFFIN_MODELS_SM_PROVIDER_H

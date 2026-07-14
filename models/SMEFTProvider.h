#ifndef GRIFFIN_MODELS_SMEFT_PROVIDER_H
#define GRIFFIN_MODELS_SMEFT_PROVIDER_H

#include "prediction/TheoryProvider.h"
#include "smeft/SMEFTInput.h"

namespace griffin {
namespace prediction {

class SMEFTProvider : public TheoryProvider {
private:
  CalculationContext context_;
  const smeft::Input* eftInput_;

public:
  SMEFTProvider(const CalculationContext& context,
                const smeft::Input& eftInput)
    : context_(context), eftInput_(&eftInput) {}

  const CalculationContext& context() const override { return context_; }
  ProviderOutputKind outputKind() const override
  {
    return ProviderOutputKind::AdditiveCorrection;
  }
  const char* name() const override { return "linear dimension-six SMEFT"; }

  const smeft::Input& eftInput() const { return *eftInput_; }

  PerturbativeSeries<Cplx> pseudoObservableSeries(
      PseudoObservable observable, Fermion fermion) const override;

  PerturbativeSeries<Cplx> amplitudeSeries(
      const ScatteringPoint& point,
      const CurrentChannel& channel) const override;
};

} // namespace prediction
} // namespace griffin

#endif // GRIFFIN_MODELS_SMEFT_PROVIDER_H

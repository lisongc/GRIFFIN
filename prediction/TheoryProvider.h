#ifndef GRIFFIN_PREDICTION_THEORY_PROVIDER_H
#define GRIFFIN_PREDICTION_THEORY_PROVIDER_H

#include "Cplx.h"
#include "classes.h"
#include "prediction/CalculationContext.h"
#include "prediction/PerturbativeSeries.h"

namespace griffin {
namespace prediction {

enum class PseudoObservable {
  SW,
  FA,
  FV,
  DeltaR,
  WMass
};

inline const char* pseudoObservableName(PseudoObservable observable)
{
  switch(observable)
  {
    case PseudoObservable::SW: return "SW";
    case PseudoObservable::FA: return "FA";
    case PseudoObservable::FV: return "FV";
    case PseudoObservable::DeltaR: return "DeltaR";
    case PseudoObservable::WMass: return "MW";
  }
  return "unknown observable";
}

struct ScatteringPoint {
  Fermion initial;
  Fermion final;
  double s;
  double cosTheta;
  double unitConversion;

  ScatteringPoint(Fermion initialState, Fermion finalState,
                  double invariantMassSquared, double cosine,
                  double conversion = 1.0)
    : initial(initialState), final(finalState),
      s(invariantMassSquared), cosTheta(cosine),
      unitConversion(conversion) {}
};

struct CurrentChannel {
  Current initial;
  Current final;

  CurrentChannel(Current initialCurrent, Current finalCurrent)
    : initial(initialCurrent), final(finalCurrent) {}
};

enum class ProviderOutputKind {
  ReferenceTheory,
  AdditiveCorrection,
  FullModel
};

class TheoryProvider {
public:
  virtual ~TheoryProvider() {}

  virtual const CalculationContext& context() const = 0;
  virtual ProviderOutputKind outputKind() const = 0;
  virtual const char* name() const = 0;

  virtual PerturbativeSeries<Cplx> pseudoObservableSeries(
      PseudoObservable observable, Fermion fermion) const = 0;

  virtual PerturbativeSeries<Cplx> amplitudeSeries(
      const ScatteringPoint& point,
      const CurrentChannel& channel) const = 0;
};

} // namespace prediction
} // namespace griffin

#endif // GRIFFIN_PREDICTION_THEORY_PROVIDER_H

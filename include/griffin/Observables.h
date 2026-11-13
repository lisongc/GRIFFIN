#ifndef GRIFFIN_PREDICTION_OBSERVABLES_H
#define GRIFFIN_PREDICTION_OBSERVABLES_H

#include "griffin/core/classes.h"

namespace griffin {
namespace prediction {

// Names accepted by the public observable interface. SMEFTCalculator rejects
// DeltaR: its finite mass-equation displacement is not a Delta-r observable.
enum class PseudoObservable {
  SW,       // effective weak mixing angle sin^2(theta_eff), not its square root
  FA,       // axial pole form factor, not the axial coupling itself
  FV,       // vector pole form factor, not the vector coupling itself
  DeltaR,
  WMass     // W pole mass, in GeV
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
  Fermion initial;       // fermion in the incoming fermion-antifermion pair
  Fermion final;         // fermion in the outgoing fermion-antifermion pair
  double s;              // center-of-mass energy squared, in GeV^2
  double cosTheta;       // cosine of the scattering angle, not the angle
  double unitConversion; // multiply GeV^-2 by this factor (0.38937966e9 for pb)

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

} // namespace prediction
} // namespace griffin

#endif // GRIFFIN_PREDICTION_OBSERVABLES_H

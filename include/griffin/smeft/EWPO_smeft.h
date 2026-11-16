/*
  File: EWPO_smeft.h
  Purpose: Define electroweak observable selections accepted by
           SMEFTCalculator. Their evaluation is in EWPO_smeft.cc;
           scattering request types are declared in xsc_smeft.h.
  Authors: Jackson Wallace and Lisong Chen
  Last updated: 2026-11-16
*/

#ifndef GRIFFIN_SMEFT_EWPO_H
#define GRIFFIN_SMEFT_EWPO_H

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

} // namespace prediction
} // namespace griffin

#endif // GRIFFIN_SMEFT_EWPO_H
